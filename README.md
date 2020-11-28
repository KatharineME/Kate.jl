[Install Docker](https://www.docker.com/get-started).

```sh
git clone https://github.com/KatharineME/analyze_sequence

cd analyze_sequence
```

Run the analyze_sequence docker container. This command makes the analyze_sequence code accessible and editable in the container. It also runs JupyterLab on port `8888` in the container and maps it to port `10000` on the host OS. Note: 'jovyan' is the default name Jupyter uses in containers.
```sh
docker run --rm --cpus=0.9 --memory=60gb -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/analyze_sequence katharineme/analyze_sequence
```

JupyterLab is running in the container at the link below. Use the `token_id` that's printed to the terminal when the container is run.

`http://127.0.0.1:10000/?token=<token_id>`