#!/bin/sh
set -eu
cd "$(dirname "$0")"
python -m pip install -r requirements.txt
python -m pytest -q
python -m streamlit run app.py --server.address=0.0.0.0 --server.port=8501
