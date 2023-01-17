#!/usr/bin/env python3

import json

class FragmentAnalyzer:

    def __init__(self, json_file: str) -> None:

        with open(json_file, "r", encoding = "utf-8") as f:
            self.data = json.load(f)
            f.close()

    def example_analysis():
        print(len(self.data))
