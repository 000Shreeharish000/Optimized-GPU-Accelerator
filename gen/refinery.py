"""Industrial model generators (literature structures, synthetic data - NOT MRPL data).

Models
------
planning   Multi-period refinery planning LP (optionally MILP): crude purchase and
           inventory, CDU with crude-specific cut yields (fixed-yield "PIMS base"
           style, crude-segregated streams so quality specs stay linear), VDU, FCC,
           reformer, diesel hydrotreater, product blending with octane / sulfur
           specs, product inventories, min/max demand with shortage penalties.
           MILP option: CDU operating-mode binaries with min-throughput and crude
           parcel integrality (cargo sizes).
           Refs: Pinto, Joly & Moro (2000); Lee, Pinto, Grossmann & Park (1996);
           review Front. Eng. Manag. (2020).
blending   Crude blending LP/MILP: meet a CDU feed spec (sulfur, API, TAN) at
           minimum cost; contract minimum lifts via binaries.
unloading  Crude-oil unloading & tank scheduling MILP (Lee et al. 1996 structure):
           vessels, storage tanks, charging tanks, CDU, time slots; transfer
           binaries with big-M ("bigm") or tight capacity-based ("tight") linking.
uc         Unit commitment MILP (3-bin formulation, min up/down, ramping, reserve).
dispatch   Economic dispatch convex QP (quadratic generator costs, network-free).

Usage:  python gen/refinery.py planning --crudes 20 --periods 52 --out data/gen/plan_20x52.mps
"""
from __future__ import annotations

import argparse
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpswriter import INF, LinearModel  # noqa: E402

CUTS = ["LPG", "NAP", "KER", "DSL", "VGO", "RES"]


def synth_assay(rng: random.Random, k: int):
    """Synthetic crude assay: lighter crudes (high API) give more light cuts."""
    api = rng.uniform(22, 44)
    sulfur = rng.uniform(0.1, 3.2)  # wt %
    light = (api - 20) / 30  # 0..0.8
    base = {
        "LPG": 0.01 + 0.02 * light,
        "NAP": 0.10 + 0.18 * light,
        "KER": 0.10 + 0.05 * light,
        "DSL": 0.20 + 0.05 * light,
        "VGO": 0.28 - 0.06 * light,
        "RES": 0.31 - 0.24 * light,
    }
    noise = {c: base[c] * rng.uniform(0.9, 1.1) for c in CUTS}
    s = sum(noise.values())
    yields = {c: v / s for c, v in noise.items()}
    cut_sulfur = {"LPG": 0.001, "NAP": 0.05, "KER": 0.15, "DSL": 0.6, "VGO": 1.2, "RES": 2.0}
    price = 70 + 0.6 * (api - 30) - 3.0 * sulfur + rng.uniform(-2, 2)  # $/bbl
    return {
        "name": f"CR{k:02d}",
        "api": api,
        "sulfur": sulfur,
        "yield": yields,
        "cut_sulfur": {c: sulfur * cut_sulfur[c] for c in CUTS},
        "nap_octane": 50 + 0.5 * (api - 20) + rng.uniform(-3, 3),
        "price": price,
    }


def planning(crudes=10, periods=12, mip=False, seed=1, out="plan.mps"):
    rng = random.Random(seed)
    C = [synth_assay(rng, k) for k in range(crudes)]
    T = list(range(periods))
    m = LinearModel(f"REFPLAN_{crudes}x{periods}{'_MIP' if mip else ''}", maximize=True)
    cdu_cap = 150.0 * max(1, crudes // 6)  # kbbl / period
    prod_price = {"LPG": 55, "GAS": 105, "JET": 100, "DSL": 98, "FO": 60}
    demand = {"LPG": 0.04, "GAS": 0.30, "JET": 0.10, "DSL": 0.35, "FO": 0.20}
    short_pen = 400.0
    for t in T:
        season = 1 + 0.1 * math.sin(2 * math.pi * t / max(periods, 1))
        # ---- crude purchase / inventory / CDU charge ----
        cdu_in = {}
        for c in C:
            n = c["name"]
            avail = rng.uniform(0.15, 0.6) * cdu_cap
            pr = c["price"] * (1 + 0.05 * math.sin(2 * math.pi * (t + rng.random()) / 12))
            if mip:
                parcel = 20.0
                np_ = m.var(f"NPAR_{n}_{t}", 0, math.floor(avail / parcel), integer=True)
                buy = m.var(f"BUY_{n}_{t}", 0, INF, cost=-pr)
                m.row(f"PARCEL_{n}_{t}", {buy: 1, np_: -parcel}, 0, 0)
            else:
                buy = m.var(f"BUY_{n}_{t}", 0, avail, cost=-pr)
            inv = m.var(f"CINV_{n}_{t}", 0, 0.5 * cdu_cap, cost=-0.3)
            run = m.var(f"CDU_{n}_{t}", 0, INF)
            prev = {f"CINV_{n}_{t - 1}": 1} if t > 0 else {}
            init = 0.05 * cdu_cap if t == 0 else 0.0
            m.row(f"CBAL_{n}_{t}", {**prev, buy: 1, inv: -1, run: -1}, -init, -init)
            cdu_in[n] = run
        cap = {r: 1 for r in cdu_in.values()}
        if mip:
            on = m.var(f"CDUMODE_{t}", 0, 1, integer=True)
            m.row(f"CDUCAP_{t}", {**cap, on: -cdu_cap}, -INF, 0)
            m.row(f"CDUMIN_{t}", {**cap, on: -0.6 * cdu_cap}, 0, INF)
            m.cost(on, -150.0)  # fixed operating cost when running
        else:
            m.row(f"CDUCAP_{t}", cap, -INF, cdu_cap)
        # ---- cuts (crude-segregated) ----
        cut = {}
        for c in C:
            n = c["name"]
            for k in CUTS:
                v = m.var(f"{k}_{n}_{t}", 0, INF)
                m.row(f"Y{k}_{n}_{t}", {v: 1, cdu_in[n]: -c["yield"][k]}, 0, 0)
                cut[(n, k)] = v
        # ---- VDU: RES -> VGO2 (0.45) + VR (0.55) ----
        vdu_cap = 0.25 * cdu_cap
        vdu_in, vgo2, vr = {}, {}, {}
        for c in C:
            n = c["name"]
            vi = m.var(f"VDUIN_{n}_{t}", 0, INF)
            rdirect = m.var(f"RESFO_{n}_{t}", 0, INF)
            m.row(f"RESSPLIT_{n}_{t}", {cut[(n, 'RES')]: 1, vi: -1, rdirect: -1}, 0, 0)
            vdu_in[n] = vi
            vgo2[n] = m.var(f"VGO2_{n}_{t}", 0, INF)
            vr[n] = m.var(f"VR_{n}_{t}", 0, INF)
            m.row(f"VDUY1_{n}_{t}", {vgo2[n]: 1, vi: -0.45}, 0, 0)
            m.row(f"VDUY2_{n}_{t}", {vr[n]: 1, vi: -0.55}, 0, 0)
        m.row(f"VDUCAP_{t}", {v: 1 for v in vdu_in.values()}, -INF, vdu_cap)
        # ---- FCC: VGO (+VGO2) -> FCCG 0.55, LCO 0.20, FLPG 0.15, SLURRY 0.07 ----
        fcc_cap = 0.35 * cdu_cap
        fcc_in = {}
        for c in C:
            n = c["name"]
            fi = m.var(f"FCCIN_{n}_{t}", 0, INF)
            m.row(f"VGOBAL_{n}_{t}", {cut[(n, 'VGO')]: 1, vgo2[n]: 1, fi: -1}, 0, 0)
            fcc_in[n] = fi
        fccg = m.var(f"FCCG_{t}", 0, INF)
        lco = m.var(f"LCO_{t}", 0, INF)
        flpg = m.var(f"FLPG_{t}", 0, INF)
        slurry = m.var(f"SLURRY_{t}", 0, INF)
        tot = {v: 1 for v in fcc_in.values()}
        m.row(f"FCCCAP_{t}", tot, -INF, fcc_cap)
        for v, y in ((fccg, 0.55), (lco, 0.20), (flpg, 0.15), (slurry, 0.07)):
            m.row(f"FCCY_{v}", {v: 1, **{k: -y for k in tot}}, 0, 0)
        m.cost(fccg, 0)
        for fi in fcc_in.values():
            m.cost(fi, -2.5)
        # LCO sulfur depends on the crude mix: bound it by the FCC feed sulfur (linear, per crude).
        # ---- Reformer: NAP -> reformate 0.85 (RON 98) ----
        ref_cap = 0.2 * cdu_cap
        ref_in = {}
        nap_gas = {}
        for c in C:
            n = c["name"]
            ri = m.var(f"REFIN_{n}_{t}", 0, INF)
            ng = m.var(f"NAPGAS_{n}_{t}", 0, INF)
            m.row(f"NAPBAL_{n}_{t}", {cut[(n, 'NAP')]: 1, ri: -1, ng: -1}, 0, 0)
            ref_in[n] = ri
            nap_gas[n] = ng
            m.cost(ri, -3.0)
        refm = m.var(f"REFORMATE_{t}", 0, INF)
        m.row(f"REFCAP_{t}", {v: 1 for v in ref_in.values()}, -INF, ref_cap)
        m.row(f"REFY_{t}", {refm: 1, **{v: -0.85 for v in ref_in.values()}}, 0, 0)
        # ---- Diesel hydrotreater: DSL + LCO -> ULSD (sulfur 0.001) ----
        dht_cap = 0.3 * cdu_cap
        dht_in, dsl_blend = {}, {}
        for c in C:
            n = c["name"]
            di = m.var(f"DHTIN_{n}_{t}", 0, INF)
            db = m.var(f"DSLDIR_{n}_{t}", 0, INF)
            m.row(f"DSLBAL_{n}_{t}", {cut[(n, 'DSL')]: 1, di: -1, db: -1}, 0, 0)
            dht_in[n] = di
            dsl_blend[n] = db
            m.cost(di, -2.0)
        lco_dht = m.var(f"LCODHT_{t}", 0, INF)
        lco_fo = m.var(f"LCOFO_{t}", 0, INF)
        m.row(f"LCOBAL_{t}", {lco: 1, lco_dht: -1, lco_fo: -1}, 0, 0)
        ulsd = m.var(f"ULSD_{t}", 0, INF)
        m.row(f"DHTCAP_{t}", {**{v: 1 for v in dht_in.values()}, lco_dht: 1}, -INF, dht_cap)
        m.row(f"DHTY_{t}", {ulsd: 1, **{v: -0.98 for v in dht_in.values()}, lco_dht: -0.95}, 0, 0)
        # ---- Product blending with specs ----
        sell = {p: m.var(f"SELL_{p}_{t}", 0, INF, cost=prod_price[p] * season) for p in prod_price}
        pinv = {p: m.var(f"PINV_{p}_{t}", 0, 0.1 * cdu_cap, cost=-0.5) for p in prod_price}
        short = {p: m.var(f"SHORT_{p}_{t}", 0, INF, cost=-short_pen) for p in prod_price}
        produced = {}
        # Gasoline: reformate (RON 98), FCC gasoline (RON 92), straight-run naphtha (crude octane); RON >= 91
        gas = m.var(f"GASPOOL_{t}", 0, INF)
        m.row(f"GASBAL_{t}", {gas: 1, refm: -1, fccg: -1, **{v: -1 for v in nap_gas.values()}}, 0, 0)
        spec = {refm: 98 - 91, fccg: 92 - 91}
        for c in C:
            spec[nap_gas[c["name"]]] = c["nap_octane"] - 91
        m.row(f"GASRON_{t}", spec, 0, INF)
        produced["GAS"] = gas
        # Jet: kerosene with sulfur <= 0.3 wt%
        jet = m.var(f"JETPOOL_{t}", 0, INF)
        kj = {}
        for c in C:
            n = c["name"]
            kj[n] = m.var(f"KERJET_{n}_{t}", 0, INF)
            kd = m.var(f"KERDSL_{n}_{t}", 0, INF)
            m.row(f"KERBAL_{n}_{t}", {cut[(n, 'KER')]: 1, kj[n]: -1, kd: -1}, 0, 0)
            dsl_blend[n + "_K"] = kd
        m.row(f"JETBAL_{t}", {jet: 1, **{v: -1 for v in kj.values()}}, 0, 0)
        m.row(f"JETSUL_{t}", {kj[c["name"]]: c["cut_sulfur"]["KER"] - 0.3 for c in C}, -INF, 0)
        produced["JET"] = jet
        # Diesel: ULSD + straight-run gasoil + kerosene, sulfur <= 0.05 wt%
        dsl = m.var(f"DSLPOOL_{t}", 0, INF)
        m.row(f"DSLPB_{t}", {dsl: 1, ulsd: -1, **{v: -1 for v in dsl_blend.values()}}, 0, 0)
        sp = {ulsd: 0.001 - 0.05}
        for c in C:
            n = c["name"]
            sp[dsl_blend[n]] = c["cut_sulfur"]["DSL"] - 0.05
            sp[dsl_blend[n + "_K"]] = c["cut_sulfur"]["KER"] - 0.05
        m.row(f"DSLSUL_{t}", sp, -INF, 0)
        produced["DSL"] = dsl
        # LPG
        lpg = m.var(f"LPGPOOL_{t}", 0, INF)
        m.row(f"LPGBAL_{t}", {lpg: 1, flpg: -1, **{cut[(c['name'], 'LPG')]: -1 for c in C}}, 0, 0)
        produced["LPG"] = lpg
        # Fuel oil: VR + direct residue + slurry + LCO; sulfur <= 3.5
        fo = m.var(f"FOPOOL_{t}", 0, INF)
        fo_terms = {fo: 1, slurry: -1, lco_fo: -1}
        fs = {slurry: 2.0 - 3.5, lco_fo: 1.0 - 3.5}
        for c in C:
            n = c["name"]
            rd = f"RESFO_{n}_{t}"
            fo_terms[vr[n]] = -1
            fo_terms[rd] = -1
            fs[vr[n]] = c["cut_sulfur"]["RES"] * 1.3 - 3.5
            fs[rd] = c["cut_sulfur"]["RES"] - 3.5
        m.row(f"FOBAL_{t}", fo_terms, 0, 0)
        m.row(f"FOSUL_{t}", fs, -INF, 0)
        produced["FO"] = fo
        # Product inventory balances and demand
        for p in prod_price:
            prev = {f"PINV_{p}_{t - 1}": 1} if t > 0 else {}
            m.row(f"PBAL_{p}_{t}", {**prev, produced[p]: 1, sell[p]: -1, pinv[p]: -1}, 0, 0)
            dmax = demand[p] * cdu_cap * season
            m.row(f"DMAX_{p}_{t}", {sell[p]: 1}, -INF, dmax * 1.2)
            m.row(f"DMIN_{p}_{t}", {sell[p]: 1, short[p]: 1}, 0.5 * dmax, INF)
    m.write(out)
    return m


def blending(crudes=12, mip=True, seed=2, out="blend.mps"):
    rng = random.Random(seed)
    C = [synth_assay(rng, k) for k in range(crudes)]
    m = LinearModel(f"CRUDEBLEND_{crudes}{'_MIP' if mip else ''}")
    D = 1000.0
    xs = {}
    for c in C:
        n = c["name"]
        avail = rng.uniform(150, 600)
        x = m.var(f"X_{n}", 0, avail, cost=c["price"])
        xs[n] = x
        if mip:
            y = m.var(f"Y_{n}", 0, 1, integer=True)
            m.row(f"LIFTMIN_{n}", {x: 1, y: -60.0}, 0, INF)   # min lift if bought
            m.row(f"LIFTMAX_{n}", {x: 1, y: -avail}, -INF, 0)
            m.cost(y, 250.0)  # contract fixed cost
    m.row("DEMAND", {x: 1 for x in xs.values()}, D, D)
    m.row("SULFUR", {xs[c["name"]]: c["sulfur"] - 1.2 for c in C}, -INF, 0)
    m.row("APILO", {xs[c["name"]]: c["api"] - 30 for c in C}, 0, INF)
    m.row("APIHI", {xs[c["name"]]: c["api"] - 38 for c in C}, -INF, 0)
    m.row("RESMAX", {xs[c["name"]]: c["yield"]["RES"] - 0.22 for c in C}, -INF, 0)
    m.write(out)
    return m


def unloading(vessels=3, storage=3, charging=2, slots=12, formulation="tight", seed=3, out="unload.mps"):
    """Crude unloading / tank scheduling MILP (component-flow linearization, Lee et al. 1996)."""
    rng = random.Random(seed)
    m = LinearModel(f"UNLOAD_V{vessels}S{storage}B{charging}T{slots}_{formulation.upper()}")
    comps = ["A", "B"]  # crude components (e.g. light sweet / heavy sour)
    sulfur = {"A": 0.2, "B": 2.5}
    V = [f"V{v}" for v in range(vessels)]
    S = [f"S{s}" for s in range(storage)]
    B = [f"B{b}" for b in range(charging)]
    T = list(range(slots))
    vol = {v: rng.uniform(60, 100) for v in V}
    arrive = {v: rng.randint(0, max(0, slots // 2 - 1)) for v in V}
    vcomp = {v: rng.choice(comps) for v in V}
    scap = {s: 150.0 for s in S}
    bcap = {b: 120.0 for b in B}
    FV, FS = 50.0, 40.0  # max transfer rates
    cdu = 30.0           # CDU demand per slot
    bsul = {B[0]: (0.8, 1.6)}
    if len(B) > 1:
        bsul[B[1]] = (1.2, 2.2)
    for b in B[2:]:
        bsul[b] = (0.5, 2.5)
    # Unloading
    for v in V:
        tot = {}
        for t in T:
            if t < arrive[v]:
                continue
            for s in S:
                x = m.var(f"XU_{v}_{s}_{t}", 0, 1, integer=True)
                f = m.var(f"FU_{v}_{s}_{t}", 0, FV)
                M = FV if formulation == "tight" else 1000.0
                m.row(f"LU_{v}_{s}_{t}", {f: 1, x: -M}, -INF, 0)
                tot[f] = 1
        m.row(f"VOL_{v}", tot, vol[v], vol[v])
        dep = m.var(f"DEP_{v}", arrive[v], slots)  # departure slot
        for t in T:
            for s in S:
                if f"XU_{v}_{s}_{t}" in m.cols:
                    m.row(f"DEPT_{v}_{s}_{t}", {dep: 1, f"XU_{v}_{s}_{t}": -(t + 1)}, 0, INF)
        m.cost(dep, 10.0)  # demurrage per slot
    # Storage balances by component
    for si, s in enumerate(S):
        for c in comps:
            for t in T:
                inv = m.var(f"IS_{s}_{c}_{t}", 0, scap[s])
                terms = {inv: 1}
                if t > 0:
                    terms[f"IS_{s}_{c}_{t - 1}"] = -1
                for v in V:
                    nm = f"FU_{v}_{s}_{t}"
                    if nm in m.cols and vcomp[v] == c:
                        terms[nm] = -1
                for b in B:
                    fc = m.var(f"FSB_{s}_{b}_{c}_{t}", 0, FS)
                    terms[fc] = 1
                init = 40.0 if c == comps[si % 2] else 0.0
                m.row(f"SB_{s}_{c}_{t}", terms, init if t == 0 else 0.0, init if t == 0 else 0.0)
            m.cost(f"IS_{s}_{c}_{slots - 1}", 0.05)
        for t in T:
            m.row(f"SCAP_{s}_{t}", {f"IS_{s}_{c}_{t}": 1 for c in comps}, -INF, scap[s])
    # Storage -> charging transfers: binary, and a tank cannot receive and deliver in the same slot
    for s in S:
        for b in B:
            for t in T:
                y = m.var(f"YSB_{s}_{b}_{t}", 0, 1, integer=True)
                M = FS if formulation == "tight" else 1000.0
                m.row(f"LSB_{s}_{b}_{t}", {**{f"FSB_{s}_{b}_{c}_{t}": 1 for c in comps}, y: -M}, -INF, 0)
                m.cost(y, 5.0)  # changeover / setup cost
    for b in B:
        for t in T:
            m.row(f"ONEIN_{b}_{t}", {f"YSB_{s}_{b}_{t}": 1 for s in S}, -INF, 1)
    # Charging tanks: component inventory, feed CDU
    for bi, b in enumerate(B):
        for t in T:
            for c in comps:
                inv = m.var(f"IB_{b}_{c}_{t}", 0, bcap[b])
                out_ = m.var(f"FB_{b}_{c}_{t}", 0, INF)
                terms = {inv: 1, out_: 1, **{f"FSB_{s}_{b}_{c}_{t}": -1 for s in S}}
                if t > 0:
                    terms[f"IB_{b}_{c}_{t - 1}"] = -1
                init = 30.0 if c == comps[bi % 2] else 10.0
                m.row(f"BB_{b}_{c}_{t}", terms, init if t == 0 else 0.0, init if t == 0 else 0.0)
            m.row(f"BCAP_{b}_{t}", {f"IB_{b}_{c}_{t}": 1 for c in comps}, -INF, bcap[b])
            lo, hi = bsul[b]
            outs = {f"FB_{b}_{c}_{t}": 1 for c in comps}
            m.row(f"BSLO_{b}_{t}", {f"FB_{b}_{c}_{t}": sulfur[c] - lo for c in comps}, 0, INF)
            m.row(f"BSHI_{b}_{t}", {f"FB_{b}_{c}_{t}": sulfur[c] - hi for c in comps}, -INF, 0)
    for t in T:
        m.row(f"CDU_{t}", {f"FB_{b}_{c}_{t}": 1 for b in B for c in comps}, cdu, cdu)
    m.write(out)
    return m


def unit_commitment(gens=20, hours=24, seed=4, out="uc.mps"):
    rng = random.Random(seed)
    m = LinearModel(f"UC_G{gens}H{hours}")
    G = []
    for g in range(gens):
        pmax = rng.uniform(50, 400)
        G.append({
            "n": f"G{g:02d}", "pmin": pmax * rng.uniform(0.3, 0.5), "pmax": pmax,
            "mc": rng.uniform(15, 60), "noload": rng.uniform(50, 400), "start": rng.uniform(200, 3000),
            "ramp": pmax * rng.uniform(0.3, 0.8), "minup": rng.randint(1, 5), "mindn": rng.randint(1, 5),
        })
    cap = sum(g["pmax"] for g in G)
    D = [0.45 * cap * (1 + 0.35 * math.sin(math.pi * (h - 6) / 12)) for h in range(hours)]
    for g in G:
        n = g["n"]
        for h in range(hours):
            u = m.var(f"U_{n}_{h}", 0, 1, integer=True)
            v = m.var(f"V_{n}_{h}", 0, 1, integer=True)
            w = m.var(f"W_{n}_{h}", 0, 1, integer=True)
            p = m.var(f"P_{n}_{h}", 0, g["pmax"], cost=g["mc"])
            m.cost(u, g["noload"])
            m.cost(v, g["start"])
            m.row(f"PMIN_{n}_{h}", {p: 1, u: -g["pmin"]}, 0, INF)
            m.row(f"PMAX_{n}_{h}", {p: 1, u: -g["pmax"]}, -INF, 0)
            prev = {f"U_{n}_{h - 1}": -1} if h > 0 else {}
            m.row(f"LOGIC_{n}_{h}", {u: 1, **prev, v: -1, w: 1}, 0, 0)
            if h > 0:
                m.row(f"RUP_{n}_{h}", {p: 1, f"P_{n}_{h - 1}": -1}, -INF, g["ramp"] + g["pmin"])
                m.row(f"RDN_{n}_{h}", {f"P_{n}_{h - 1}": 1, p: -1}, -INF, g["ramp"] + g["pmin"])
            ups = {f"V_{n}_{k}": 1 for k in range(max(0, h - g["minup"] + 1), h + 1)}
            m.row(f"MINUP_{n}_{h}", {**ups, u: -1}, -INF, 0)
            dns = {f"W_{n}_{k}": 1 for k in range(max(0, h - g["mindn"] + 1), h + 1)}
            m.row(f"MINDN_{n}_{h}", {**dns, u: 1}, -INF, 1)
    for h in range(hours):
        m.row(f"DEM_{h}", {f"P_{g['n']}_{h}": 1 for g in G}, D[h], INF)
        m.row(f"RES_{h}", {f"U_{g['n']}_{h}": g["pmax"] for g in G}, 1.1 * D[h], INF)
    m.write(out)
    return m


def dispatch(gens=50, hours=24, seed=5, out="dispatch.qps"):
    rng = random.Random(seed)
    m = LinearModel(f"EDQP_G{gens}H{hours}")
    G = [{"n": f"G{g:02d}", "pmax": rng.uniform(50, 400), "a": rng.uniform(10, 40), "b": rng.uniform(0.002, 0.05),
          "ramp": rng.uniform(40, 150)} for g in range(gens)]
    cap = sum(g["pmax"] for g in G)
    for h in range(hours):
        d = 0.55 * cap * (1 + 0.3 * math.sin(math.pi * (h - 6) / 12))
        for g in G:
            p = m.var(f"P_{g['n']}_{h}", 0, g["pmax"], cost=g["a"])
            m.qterm(p, p, 2 * g["b"])  # 1/2 * (2b) p^2 = b p^2
            if h > 0:
                m.row(f"R_{g['n']}_{h}", {p: 1, f"P_{g['n']}_{h - 1}": -1}, -g["ramp"], g["ramp"])
        m.row(f"D_{h}", {f"P_{g['n']}_{h}": 1 for g in G}, d, d)
    m.write(out)
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model", choices=["planning", "blending", "unloading", "uc", "dispatch", "suite"])
    ap.add_argument("--crudes", type=int, default=10)
    ap.add_argument("--periods", type=int, default=12)
    ap.add_argument("--mip", action="store_true")
    ap.add_argument("--formulation", default="tight", choices=["tight", "bigm"])
    ap.add_argument("--gens", type=int, default=20)
    ap.add_argument("--hours", type=int, default=24)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.model == "suite":
        d = a.out or "data/gen"
        os.makedirs(d, exist_ok=True)
        for c, t in [(6, 4), (10, 12), (20, 12), (20, 52), (30, 52), (40, 104), (60, 156)]:
            m = planning(c, t, False, a.seed, f"{d}/plan_{c}x{t}.mps")
            print(f"{m.name}: {m.num_rows} rows, {m.num_cols} cols")
        for c, t in [(6, 4), (10, 6), (10, 12)]:
            m = planning(c, t, True, a.seed, f"{d}/plan_{c}x{t}_mip.mps")
            print(f"{m.name}: {m.num_rows} rows, {m.num_cols} cols")
        for f in ("tight", "bigm"):
            m = unloading(3, 3, 2, 12, f, a.seed, f"{d}/unload_{f}.mps")
            print(f"{m.name}: {m.num_rows} rows, {m.num_cols} cols")
        for k, mip in ((12, False), (12, True), (30, True)):
            m = blending(k, mip, a.seed, f"{d}/blend_{k}{'_mip' if mip else ''}.mps")
            print(f"{m.name}: {m.num_rows} rows, {m.num_cols} cols")
        for g, h in ((10, 24), (30, 24)):
            m = unit_commitment(g, h, a.seed, f"{d}/uc_{g}x{h}.mps")
            print(f"{m.name}: {m.num_rows} rows, {m.num_cols} cols")
        m = dispatch(50, 24, a.seed, f"{d}/dispatch_50x24.qps")
        print(f"{m.name}: {m.num_rows} rows, {m.num_cols} cols")
        return
    out = a.out or f"{a.model}.mps"
    if a.model == "planning":
        m = planning(a.crudes, a.periods, a.mip, a.seed, out)
    elif a.model == "blending":
        m = blending(a.crudes, a.mip, a.seed, out)
    elif a.model == "unloading":
        m = unloading(formulation=a.formulation, seed=a.seed, out=out)
    elif a.model == "uc":
        m = unit_commitment(a.gens, a.hours, a.seed, out)
    else:
        m = dispatch(a.gens, a.hours, a.seed, out)
    print(f"{m.name}: {m.num_rows} rows, {m.num_cols} cols -> {out}")


if __name__ == "__main__":
    main()
