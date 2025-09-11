# 🔬 Sanger Assembly App

A lightweight, Python-based Streamlit web application for generating high-quality consensus sequences from paired-end Sanger `.ab1` chromatogram files.

Designed for use in molecular biology, diagnostics, and teaching labs.

---

## ✨ Features

- 🧬 Input paired `.ab1` Sanger reads (forward and reverse)
  <img width="1225" height="603" alt="image" src="https://github.com/user-attachments/assets/b0f24f27-62cc-409e-a171-4df269c2c16b" />
- ✂️ Automatic quality-based trimming\
  Analysis Parameters Default:\
  Quality Cutoff: 20 (Minimum phred score threshold)\
  Window Size: 10 (Size of sliding window for quality assessment)\
  Score Margin: 1.5 (How much better must the Best Window be to override the First Acceptable Window? If the score margin is 0, select the First Acceptable Window.)
  
  <img width="475" height="574" alt="image" src="https://github.com/user-attachments/assets/727d5460-36fa-44b8-aa01-5240f01a7d39" />
- 🔁 Reverse complement + reverse Phred score handling
  <img width="1202" height="545" alt="image" src="https://github.com/user-attachments/assets/35dfe479-0a69-47d1-89c2-8d9266e9af70" />
- 🧩 Overlap alignment
  <img width="1384" height="772" alt="image" src="https://github.com/user-attachments/assets/67a26a61-6645-4c57-bb95-c121f34e248d" />
- 📊 Output stats
  Metric  Value
  Matches  277
  Mismatches  12
  Gaps  983
  Alignment Score  227.5
  Consensus_Length  1272
  Quality_Cutoff  20
  Window_Size  10
  Score_Margin  1.5

- 💾 Export to FASTA, CSV, and optional BLAST-ready format
- 🌐 Web interface (Streamlit) — runs locally or on Streamlit Cloud

---

## 🚀 Getting Started

### 📦 Install requirements

```bash
pip install -r requirements.txt
