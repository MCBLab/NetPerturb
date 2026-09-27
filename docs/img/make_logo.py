"""Build the NetPerturb hex logo as SVG (text converted to paths) plus PNGs.

The title is drawn from Atkinson Hyperlegible glyph outlines rather than as
<text>, so the SVG renders identically in browsers, rasterisers and anywhere
the font is not installed. Regenerate from the repository root with:

    pip install fonttools brotli cairosvg pillow
    python docs/img/make_logo.py \
        docs/fonts/Atkinson-Hyperlegible-Bold-102a.woff2 docs/img \
        docs/fonts/Atkinson-Hyperlegible-Regular-102a.woff2

Colours are the MCB Lab site's Sandstone palette, plus a red for the
perturbed node.
"""
import math, sys
from fontTools.ttLib import TTFont
from fontTools.pens.svgPathPen import SVGPathPen
from fontTools.pens.transformPen import TransformPen
import cairosvg

FONT = sys.argv[1]
OUT = sys.argv[2]  # output directory

BLUE, CREAM, DARK, SAND, RED = "#325d88", "#f8f5f0", "#3e3f3a", "#dfd7ca", "#c8352b"


def text_path(text, size, cx, baseline, font):
    gs = font.getGlyphSet()
    cmap = font.getBestCmap()
    upm = font["head"].unitsPerEm
    scale = size / upm
    names = [cmap[ord(c)] for c in text]
    width = sum(gs[n].width for n in names) * scale
    x = cx - width / 2
    d = []
    for n in names:
        pen = SVGPathPen(gs)
        gs[n].draw(TransformPen(pen, (scale, 0, 0, -scale, x, baseline)))
        d.append(pen.getCommands())
        x += gs[n].width * scale
    return " ".join(d)


def hexagon(cx, cy, r):
    pts = [(cx + r * math.cos(math.radians(a)), cy + r * math.sin(math.radians(a)))
           for a in range(-90, 270, 60)]
    return " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)


font = TTFont(FONT)
W, H = 400, 462
cx, cy = W / 2, H / 2

# network: the perturbed (red) node, its partners, and the edges among them
red = (184, 168)
nodes = {"a": (112, 112), "b": (226, 74), "c": (296, 150), "d": (104, 220), "e": (206, 244)}
solid = [("a", "b"), ("b", "c"), ("d", "e"), ("a", "d")]

edges = "\n".join(
    f'<line x1="{nodes[p][0]}" y1="{nodes[p][1]}" x2="{nodes[q][0]}" y2="{nodes[q][1]}" '
    f'stroke="{DARK}" stroke-width="6" stroke-linecap="round"/>' for p, q in solid)
dashed = "\n".join(
    f'<line x1="{red[0]}" y1="{red[1]}" x2="{x}" y2="{y}" stroke="{RED}" stroke-width="5" '
    f'stroke-dasharray="11 9" stroke-linecap="round"/>' for x, y in nodes.values())
dots = "\n".join(
    f'<circle cx="{x}" cy="{y}" r="15" fill="{DARK}" stroke="{CREAM}" stroke-width="4"/>'
    for x, y in nodes.values())

# a small cell to the lower right: membrane, nucleus, a few organelles
cell_x, cell_y = 300, 236
cell = f'''
<g>
  <ellipse cx="{cell_x}" cy="{cell_y}" rx="38" ry="34" fill="#e3ebf3" stroke="{BLUE}" stroke-width="5"/>
  <circle cx="{cell_x + 6}" cy="{cell_y - 4}" r="14" fill="{BLUE}"/>
  <circle cx="{cell_x - 22}" cy="{cell_y + 14}" r="4.5" fill="{BLUE}" opacity="0.55"/>
  <circle cx="{cell_x + 26}" cy="{cell_y + 20}" r="3.5" fill="{BLUE}" opacity="0.55"/>
  <circle cx="{cell_x - 18}" cy="{cell_y - 20}" r="3" fill="{BLUE}" opacity="0.55"/>
</g>'''

title = text_path("NetPerturb", 54, cx, 330, font)
REG = TTFont(sys.argv[3])
tag = text_path("single-cell perturbation", 17, cx, 384, REG)

svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" role="img" aria-label="NetPerturb">
<title>NetPerturb</title>
<polygon points="{hexagon(cx, cy, 222)}" fill="{BLUE}"/>
<polygon points="{hexagon(cx, cy, 206)}" fill="{CREAM}"/>
{edges}
{dashed}
<circle cx="{red[0]}" cy="{red[1]}" r="34" fill="none" stroke="{RED}" stroke-width="3" opacity="0.35"/>
{dots}
<circle cx="{red[0]}" cy="{red[1]}" r="22" fill="{RED}" stroke="{CREAM}" stroke-width="4"/>
{cell}
<path d="{title}" fill="{DARK}"/>
<rect x="140" y="346" width="120" height="6" rx="3" fill="{RED}"/>
<path d="{tag}" fill="#5f5d58"/>
</svg>
'''

open(f"{OUT}/logo.svg", "w").write(svg)
cairosvg.svg2png(bytestring=svg.encode(), write_to=f"{OUT}/logo.png", output_width=600)

# Text-free mark for favicons: the network and cell recentred in the hexagon,
# since the title is unreadable at 16-32 px.
mark = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">
<polygon points="{hexagon(cx, cy, 222)}" fill="{BLUE}"/>
<polygon points="{hexagon(cx, cy, 200)}" fill="{CREAM}"/>
<g transform="translate(200 231) scale(1.12) translate(-203 -160)">
{edges}
{dashed}
{dots}
<circle cx="{red[0]}" cy="{red[1]}" r="24" fill="{RED}" stroke="{CREAM}" stroke-width="4"/>
{cell}
</g>
</svg>
'''
open(f"{OUT}/mark.svg", "w").write(mark)
from PIL import Image
import io
for name, size in [("favicon-16x16.png", 16), ("favicon-32x32.png", 32), ("apple-touch-icon.png", 180),
                   ("android-chrome-192x192.png", 192), ("android-chrome-512x512.png", 512)]:
    png = cairosvg.svg2png(bytestring=mark.encode(), output_width=round(size * W / H), output_height=size)
    img = Image.open(io.BytesIO(png)).convert("RGBA")
    canvas = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    canvas.paste(img, ((size - img.width) // 2, 0), img)
    canvas.save(f"{OUT}/{name}")
Image.open(f"{OUT}/android-chrome-512x512.png").save(f"{OUT}/favicon.ico", sizes=[(16, 16), (32, 32), (48, 48)])
