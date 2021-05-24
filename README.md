[Install Docker](https://www.docker.com/get-started).

```sh
git clone https://github.com/KatharineME/ProcessSequence.jl

cd ProcessSequence.jl
```

Run the [process_sequence](https://hub.docker.com/repository/docker/katharineme/process_sequence) docker container. This command makes the `ProcessSequence.jl` code accessible and editable in the container. It also runs JupyterLab on port `8888` in the container and maps it to port `10000` on the host OS. Note: 'jovyan' is the default name Jupyter uses in containers.
```sh
docker run --rm -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/ProcessSequence.jl katharineme/process_sequence
```

JupyterLab is now running at the link below. Use the `token_id` that was printed to the terminal.

`http://127.0.0.1:10000/?token=<token_id>`
