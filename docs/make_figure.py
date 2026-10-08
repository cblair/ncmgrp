"""Generate docs/space_use.png from the collar data in this repository.

Estimates the Brownian motion variance by maximum likelihood, builds the
Brownian bridge utilisation distribution over a raster grid, and renders it
beside the raw fixes and the example covariate grid.

Run from anywhere:

    python3 docs/make_figure.py

Requires numpy and Pillow. This is a documentation figure only. It is not part
of the analysis pipeline, which is the R code under BB/ and SYN2/.
"""
import math, os
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIX  = os.path.join(ROOT, "BB", "data", "003SCF.dat")
GRID = os.path.join(ROOT, "BB", "data", "003SCF.datexample_grid.asc")
OUT  = os.path.join(ROOT, "docs", "space_use.png")
FONTDIR = "/System/Library/Fonts/Supplemental/"   # adjust on non-macOS

DELTA2 = 20.0**2     # GPS collar location error, metres, assumed
MAXGAP = 1500.0      # minutes; segments longer than this are collar gaps, dropped
PX     = 560         # raster side, pixels

# ---------------------------------------------------------------- load fixes
rows = [l.split("\t") for l in open(FIX).read().splitlines() if l.strip()][1:]
t, X, Y = [], [], []
for r in rows:
    mm, dd, yy = (int(v) for v in r[0].split("/"))
    hh, mi = (int(v) for v in r[1].split(":"))
    a = (14 - mm)//12; y2 = yy + 4800 - a; m2 = mm + 12*a - 3
    jdn = dd + (153*m2+2)//5 + 365*y2 + y2//4 - y2//100 + y2//400 - 32045
    t.append(jdn*1440.0 + hh*60 + mi); X.append(float(r[2])); Y.append(float(r[3]))
o = np.argsort(t)
t, X, Y = np.array(t)[o], np.array(X)[o], np.array(Y)[o]
print(f"fixes={len(t)}  span={(t[-1]-t[0])/1440:.0f} d  median interval={np.median(np.diff(t)):.0f} min")

# ------------------------------- ML Brownian motion variance (Horne et al. 2007)
def negll(s2m):
    """Treat each interior fix as unobserved; score it against the bridge
    spanning its two neighbours."""
    T = t[2:] - t[:-2]
    ok = (T > 0) & (T < 2*MAXGAP); T = T[ok]
    al = ((t[1:-1]-t[:-2])[ok]) / T
    mux = X[:-2][ok] + al*(X[2:][ok]-X[:-2][ok])
    muy = Y[:-2][ok] + al*(Y[2:][ok]-Y[:-2][ok])
    v   = T*al*(1-al)*s2m + ((1-al)**2 + al**2)*DELTA2
    d2  = (X[1:-1][ok]-mux)**2 + (Y[1:-1][ok]-muy)**2
    return np.sum(np.log(2*math.pi*v) + d2/(2*v))

gr = (math.sqrt(5)-1)/2
a_, b_ = math.log(1e-4), math.log(1e4)
c_, d_ = b_-gr*(b_-a_), a_+gr*(b_-a_)
for _ in range(90):
    if negll(math.exp(c_)) < negll(math.exp(d_)): b_, d_, c_ = d_, c_, b_-gr*(d_-a_)
    else:                                         a_, c_, d_ = c_, d_, a_+gr*(b_-c_)
S2M = math.exp((a_+b_)/2)
print(f"ML sigma_m^2 = {S2M:.2f} m^2/min")

# ------------------------------------------------- Brownian bridge UD on a grid
pad = 500.0
x0, x1, y0, y1 = X.min()-pad, X.max()+pad, Y.min()-pad, Y.max()+pad
span = max(x1-x0, y1-y0)
cx, cy = (x0+x1)/2, (y0+y1)/2
x0, y1 = cx-span/2, cy+span/2
res = span/PX
gx = x0 + (np.arange(PX)+0.5)*res
gy = y1 - (np.arange(PX)+0.5)*res          # row 0 is north
UD = np.zeros((PX, PX))

for i in range(len(t)-1):
    T = t[i+1]-t[i]
    if T <= 0 or T > MAXGAP: continue
    # step count scales with bridge length, so long segments smear smoothly
    # instead of beading into visibly discrete kernels
    NA = int(min(400, max(24, math.hypot(X[i+1]-X[i], Y[i+1]-Y[i])/res*2.5)))
    for k in range(NA):
        al = (k+0.5)/NA
        v  = T*al*(1-al)*S2M + ((1-al)**2 + al**2)*DELTA2
        mx = X[i] + al*(X[i+1]-X[i]); my = Y[i] + al*(Y[i+1]-Y[i])
        r  = 3.5*math.sqrt(v)
        ia = max(0, int((mx-r-x0)/res)); ib = min(PX, int((mx+r-x0)/res)+1)
        ja = max(0, int((y1-(my+r))/res)); jb = min(PX, int((y1-(my-r))/res)+1)
        if ib <= ia or jb <= ja: continue
        dx = gx[ia:ib]-mx; dy = gy[ja:jb]-my
        UD[ja:jb, ia:ib] += (T/NA)/(2*math.pi*v)*np.exp(-(dy[:,None]**2+dx[None,:]**2)/(2*v))
UD /= UD.sum()

def vol_mask(p, pct):
    f = np.sort(p.ravel())[::-1]; c = np.cumsum(f)
    return p >= f[np.searchsorted(c, pct*c[-1])]

def outline(m):
    e = m.copy()
    for s, ax in ((1,0),(-1,0),(1,1),(-1,1)): e &= np.roll(m, s, axis=ax)
    return m & ~e

c50, c95 = outline(vol_mask(UD,0.50)), outline(vol_mask(UD,0.95))
A50 = vol_mask(UD,0.50).sum()*res*res/1e6
A95 = vol_mask(UD,0.95).sum()*res*res/1e6
print(f"50% UD {A50:.2f} km2   95% UD {A95:.2f} km2")

# ------------------------------------------------------ covariate grid, resampled
ls = open(GRID).read().splitlines()
hd = {ls[i].split()[0].upper(): float(ls[i].split()[1]) for i in range(6)}
G  = np.array([[float(v) for v in l.split()] for l in ls[6:] if l.strip()])
nc, nr = int(hd["NCOLS"]), int(hd["NROWS"])
gx0, gy0, cs = hd["XLLCORNER"], hd["YLLCORNER"], hd["CELLSIZE"]
cov = np.full((PX,PX), np.nan)
for j in range(PX):
    rr = int((gy0 + nr*cs - (y1-(j+0.5)*res))/cs)
    if not (0 <= rr < nr): continue
    cc = ((x0 + (np.arange(PX)+0.5)*res - gx0)/cs).astype(int)
    ok = (cc >= 0) & (cc < nc)
    cov[j, ok] = G[rr, cc[ok]]
cov = np.where(cov <= 0, np.nan, cov)

# ------------------------------------------------------------------- rendering
BG, BORDER, INK, DIM = (10,10,10), (42,42,42), (236,236,236), (132,132,132)
f_t  = ImageFont.truetype(FONTDIR+"Arial Bold.ttf", 23)
f_h  = ImageFont.truetype(FONTDIR+"Arial Bold.ttf", 17)
f_s  = ImageFont.truetype(FONTDIR+"Arial.ttf", 14)
f_xs = ImageFont.truetype(FONTDIR+"Arial.ttf", 12)

def to_px(x, y): return (x-x0)/res, (y1-y)/res
def ramp(a, lo, hi, g): return np.clip((a-lo)/(hi-lo+1e-30), 0, 1)**g

def panel(kind):
    if kind == "fix":
        im = Image.new("RGB", (PX,PX), BG); d = ImageDraw.Draw(im, "RGBA")
        for i in range(len(X)-1):
            d.line([to_px(X[i],Y[i]), to_px(X[i+1],Y[i+1])], fill=(150,150,150,26))
        for i in range(len(X)):
            a,b = to_px(X[i],Y[i]); d.ellipse([a-1.5,b-1.5,a+1.5,b+1.5], fill=(248,248,248,170))
        return im
    if kind == "ud":
        v = ramp(UD, 0, np.percentile(UD,99.92), 0.42)
        img = np.dstack([v,v,v])*0.97 + np.array(BG)/255.0*(1-v[...,None])
        im = Image.fromarray((np.clip(img,0,1)*255).astype(np.uint8)); d = ImageDraw.Draw(im,"RGBA")
        for m, al in ((c95,120),(c50,235)):
            ys,xs = np.where(m); d.point(list(zip(xs.tolist(),ys.tolist())), fill=(255,255,255,al))
        return im
    m = ~np.isnan(cov); v = np.zeros((PX,PX))
    if m.any():
        v[m] = ramp(cov[m], np.nanpercentile(cov,2), np.nanpercentile(cov,98), 0.85)*0.52
    img = np.dstack([v,v,v]) + np.array(BG)/255.0*(1-v[...,None])
    im = Image.fromarray((np.clip(img,0,1)*255).astype(np.uint8)); d = ImageDraw.Draw(im,"RGBA")
    for lv, al in ((0.30,70),(0.60,100),(0.85,130)):
        mm = v >= lv*0.52
        e = mm.copy()
        for s,ax in ((1,0),(-1,0),(1,1),(-1,1)): e &= np.roll(mm,s,axis=ax)
        ys,xs = np.where(mm & ~e); d.point(list(zip(xs.tolist(),ys.tolist())), fill=(255,255,255,al))
    for mk, al in ((c95,200),(c50,255)):
        ys,xs = np.where(mk); d.point(list(zip(xs.tolist(),ys.tolist())), fill=(255,255,255,al))
    return im

GAP, ML, TOP = 44, 54, 92
W, H = ML*2 + PX*3 + GAP*2, TOP + PX + 128
sheet = Image.new("RGB", (W,H), BG); dr = ImageDraw.Draw(sheet)
dr.text((ML,26), "Brownian bridge space use, collar 003SCF", font=f_t, fill=INK)
dr.text((ML,58), f"{len(t)} GPS fixes over {(t[-1]-t[0])/1440:.0f} days at 5 h intervals"
                 "  ·  North Cascade Mountain Goat Research Project", font=f_s, fill=DIM)

caps = [("Observed fixes",      ["Where the collar reported.", "Gaps between fixes are unobserved."]),
        ("Brownian bridge UD",  ["Movement-model probability surface.", "Contours: 50% and 95% volume."]),
        ("Habitat covariate",   ["Example grid shipped with the repo,", "30 m cells, with UD contours over it."])]
for i, (ttl, sub) in enumerate(caps):
    x = ML + i*(PX+GAP)
    sheet.paste(panel(["fix","ud","cov"][i]), (x,TOP))
    dr.rectangle([x,TOP,x+PX-1,TOP+PX-1], outline=BORDER)
    dr.text((x,TOP+PX+16), ttl, font=f_h, fill=INK)
    for k,s in enumerate(sub): dr.text((x,TOP+PX+42+k*19), s, font=f_xs, fill=DIM)

bx, by, bl = ML+16, TOP+PX-22, 500.0/res
dr.line([(bx,by),(bx+bl,by)], fill=INK, width=2)
dr.text((bx,by-20), "500 m", font=f_xs, fill=INK)
dr.text((ML,H-26), f"sigma_m^2 = {S2M:.2f} m^2/min by maximum likelihood  ·  50% UD {A50:.2f} km^2"
                   f"  ·  95% UD {A95:.2f} km^2  ·  assumed collar error 20 m", font=f_xs, fill=(104,104,104))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
sheet.save(OUT)
print("wrote", OUT, sheet.size)
