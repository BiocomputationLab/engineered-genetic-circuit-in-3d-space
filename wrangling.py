import openpyxl
import pathlib

libraries = ["4T5T", "5L14L3", "6L3", "TS2MmeI", "5A4A", "6A", "6T", "TS2dm0"]

def library_data(lib):
    path = pathlib.Path(__file__)
    return openpyxl.load_workbook(filename = f"{path.parent}/data/KT2440_{lib}.xlsx")
