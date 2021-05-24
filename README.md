#### 1. [Install Docker](https://www.docker.com/get-started)

#### 2. Get ProcessSequence.jl

```sh
git clone https://github.com/KatharineME/ProcessSequence.jl

cd ProcessSequence.jl
```

#### 3. Run the [process_sequence](https://hub.docker.com/repository/docker/katharineme/process_sequence) docker container

```sh
docker run --rm -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/ProcessSequence.jl katharineme/process_sequence
```

This command makes the `ProcessSequence.jl` code accessible and editable in the container. It also runs JupyterLab on port `8888` in the container and maps it to port `10000` on the host OS. Note: 'jovyan' is the default name Jupyter uses in containers. 

#### 4. Access JupyterLab at this link:

`http://127.0.0.1:10000/?token=<token_id>`

Use the `token_id` that was printed to the terminal.
