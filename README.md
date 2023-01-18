# eubic23_internal_ions

## Usage json_to_dataframes

```python
import pandas as pd
from json_to_dataframes import json_to_dataframes

dataframes = json_to_dataframes("data.json")

fragments_dataframe = dataframes[0]
spectra_dataframe = dataframes[1]

fragments_dataframe.head()
```
