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

## Usage web app

Maximum upload file size will be 5GB
```
streamlit run streamlit_app.py --server.maxUploadSize 5000
```

## Quick Setup with Docker

- Install [Docker](https://docs.docker.com/engine/install/).
- To run Fragment Explorer on your own server:
  ```bash
  docker run -d --restart always -p 80:8501 michabirklbauer/fragmentexplorer:latest
  ```
- To run Fragment Explorer locally:
  ```bash
  docker run -p 8501:8501 michabirklbauer/fragmentexplorer:latest
  ```
