#!/usr/bin/env python3

import json
from typing import Union
from typing import TextIO

class FragmentAnalyzer:

    def __init__(self, json_file: Union[str, TextIO], file_object: bool = False) -> None:

            if file_object:
                self.data = json.load(json_file)
            else:
                with open(json_file, "r", encoding = "utf-8") as f:
                    self.data = json.load(f)
                    f.close()

    def example_analysis(self):
        print(len(self.data))

def main():

    fa = FragmentAnalyzer("sample.json")
    fa.example_analysis()

if __name__ == "__main__":
    main()
