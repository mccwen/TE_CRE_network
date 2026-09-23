from pathlib import Path
import random

import pandas as pd
from fuc import pybed
from pybedtools import BedTool

from TE_CRE_network import steamer_enhancer as se


DATA_DIR = Path(__file__).resolve().parent / "test_data"


def test_get_enhancers_parses_cicero_peak_pairs():
    result = se.get_enhancers(DATA_DIR / "test_sig_enhancer.csv")

    assert list(result.columns) == ["chr", "start_position", "end_position"]
    assert len(result) == 10
    assert set(result["chr"]) == {"chr2"}
    assert pd.api.types.is_integer_dtype(result["start_position"])
    assert pd.api.types.is_integer_dtype(result["end_position"])
    assert (result["start_position"] < result["end_position"]).all()


def test_create_bed_for_tes_normalizes_coordinates(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    result = se.create_bed_for_TEs(DATA_DIR / "test.tsv")

    expected = pd.read_csv(DATA_DIR / "test.bed", sep="\t")
    expected_bed = pybed.BedFrame.from_frame(meta=[], data=expected)

    assert result.to_string() == expected_bed.to_string()
    assert (tmp_path / "TEs.bed").is_file()


def test_get_nearby_enhancers_filters_chromosome_and_window(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    random.seed(2017)
    enhancers = pd.DataFrame(
        {
            "chr": ["chr2", "chr2", "chr1"],
            "start_position": [100, 2_000_000, 100],
            "end_position": [200, 2_000_100, 200],
        }
    )

    _, barcodes = se.get_nearby_enhancers(enhancers, "chr2", 150, 175)

    assert len(barcodes) == 1
    assert barcodes.str.fullmatch(r"[A-Za-z0-9]{10}").all()
    fragment_path = tmp_path / "Frag.bed"
    assert fragment_path.is_file()
    fields = fragment_path.read_text().strip().split("\t")
    assert fields[:3] == ["chr2", "100", "200"]
    assert fields[3] == barcodes.iloc[0]


def test_make_cell_x_element_matrix_uses_intersection_data():
    intersection_path = DATA_DIR / "test_BEDintersection_data.bed"
    expected = pd.read_csv(DATA_DIR / "test_TE_barcode.csv")

    unique_table, family_table, unique_index, family_index, barcode_index = (
        se.make_cell_x_element_matrix(
            BedTool(str(intersection_path)),
            expected["barcode"],
        )
    )

    assert list(family_index) == expected["TE"].tolist()
    assert list(barcode_index) == expected["barcode"].tolist()
    assert len(unique_index) == len(expected)
    assert unique_table["data"].sum() == len(expected)
    assert family_table["data"].sum() == len(expected)
    assert len(unique_table) == len(expected)
    assert len(family_table) == len(expected)
