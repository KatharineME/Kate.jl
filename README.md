#### 1. [Install Docker](https://www.docker.com/get-started)

#### 2. Get Kate.jl

```sh
git clone https://github.com/KatharineME/Kate.jl

cd Kate.jl
```

#### 3. Run the [kate](https://hub.docker.com/repository/docker/katharineme/kate) Docker container

```sh
docker run --rm -p 10000:8888 -e JUPYTER_ENABLE_LAB=yes -v "$PWD":/home/jovyan/Kate.jl katharineme/kate
```

This command makes the `Kate.jl` code accessible and editable in the container. It also runs JupyterLab on port `8888` in the container and maps it to port `10000` on the host OS. Note: 'jovyan' is the default name Jupyter uses in containers. 

Learn more about Docker.

#### 4. Access JupyterLab at this link

`http://127.0.0.1:10000/?token=<token_id>`

Use the `token_id` that was printed to the terminal.

#### Optional: Run commands in container's terminal

```sh
docker exec -it <container_id> /bin/bash
```
