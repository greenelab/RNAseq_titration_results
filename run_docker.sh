# script to get docker image up and running

sudo snap install docker

sleep 3 && sudo chmod 666 /var/run/docker.sock

docker pull envest/rnaseq_titration_results:R-4.1.2

docker run -d --rm --mount type=volume,dst=/home/rstudio,volume-driver=local,volume-opt=type=none,volume-opt=o=bind,volume-opt=device=/home/ubuntu -e PASSWORD=onecupatatime envest/rnaseq_titration_results:R-4.1.2

container_id=$(docker ps -ql)

docker exec -it $container_id bash
