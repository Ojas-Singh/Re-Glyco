# Information



# Installation

```
sudo apt install build-essential
conda create -n reglyco python=3.12
conda activate reglyco
pip install -r requirements.txt

curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustup update

cd glycors
maturin develop --release
```

python3 -m pip -r requirements.txt
export PATH="$HOME/.local/bin:$PATH"


modify config.py 
run

python main.py

gunicorn -w 4 api.py:app

flask -A api.py run     


# things for oracle for Hosting Websites
```
sudo iptables -P INPUT ACCEPT
sudo iptables -P OUTPUT ACCEPT
sudo iptables -P FORWARD ACCEPT
sudo iptables -F
sudo ufw allow 22
sudo ufw allow 80
sudo ufw allow 8080
sudo ufw allow 443
sudo ufw enable
```

localhost/loopback
```
sudo iptables -t nat -I OUTPUT -p tcp -d 127.0.0.1 --dport 80 -j REDIRECT --to-ports 3000
```
external
```
sudo iptables -t nat -I PREROUTING -p tcp --dport 80 -j REDIRECT --to-ports 3000
```



