## My personal repo for Julia Docker

*Steps*

* To build a container ```bash compose up -d ```

* To build an image ```bash  docker build -t mmogib/myjulia:latest --build-arg JULIA_VERSION=1.9.0 . ```

* Run ```docker exec -it [name] zsh```


