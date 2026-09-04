"""Read-only PDF-content checks and generated contact sheets for visual QA."""

from pathlib import Path
import subprocess

from PIL import Image, ImageDraw
from pypdf import PdfReader


ROOT = Path(__file__).resolve().parent
PDF = ROOT / "output/pdf/explicit_chen_c249_final.pdf"
TMP = ROOT / "tmp/pdfs"
reader = PdfReader(PDF)
text = "\n".join(page.extract_text() for page in reader.pages)
if "??" in text:
    raise ValueError("An unresolved reference appears in the PDF")
print("PDF pages:", len(reader.pages))
for i, page in enumerate(reader.pages, 1):
    content = page.extract_text()
    if len(content.strip()) < 50:
        raise ValueError(f"Unexpected empty/sparse page {i}")
    if "external-input condition ledger" in content.lower():
        print("External-input ledger page:", i)
    if "one forward ledger" in content.lower():
        print("Forward-ledger page:", i)
TMP.mkdir(parents=True, exist_ok=True)
subprocess.run(
    ["pdftoppm", "-r", "100", "-png", str(PDF), str(TMP / "c249")],
    check=True,
)
pages = [TMP / f"c249-{i:02d}.png" for i in range(1, len(reader.pages) + 1)]
for group in range((len(pages) + 5) // 6):
    sheet = Image.new("RGB", (1470, 1470), "#e6e8eb")
    draw = ImageDraw.Draw(sheet)
    for j, path in enumerate(pages[group * 6:(group + 1) * 6]):
        with Image.open(path) as source:
            source.thumbnail((476, 696))
            x, y = (j % 3) * 490 + 7, (j // 3) * 735 + 28
            sheet.paste(source, (x, y))
            draw.text((x, y - 18), f"Page {group * 6 + j + 1}", fill="black")
    output = TMP / f"c249-contact-{group + 1}.png"
    sheet.save(output)
    print(output)

