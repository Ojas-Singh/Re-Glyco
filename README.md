# Important!
ASN-NAG

From Privateer 
Target Phi = -97.5, target Psi = 178

Need Verification.
Phi = O5 C1 ND2 CG
        CG-ND2-C1-O5

Psi = C1 ND2 CG CB
        CB-CG-ND2-C1


To-do
1. Clean code with phi psi flex
2. Write a test for AlphaFold
3. Write a flip attach


# Installation


```


conda create -n reglyco python=3.10
conda activate reglyco
pip install -r requirements.txt

sudo apt install build-essential

cd glycors
maturin develop --release
cd ..
streamlit run main.py
```

python3 -m pip -r requirements.txt
export PATH="$HOME/.local/bin:$PATH"


modify config.py 
run

python main.py






# things for oracle for Hosting Websites
sudo iptables -P INPUT ACCEPT
sudo iptables -P OUTPUT ACCEPT
sudo iptables -P FORWARD ACCEPT
sudo iptables -F
sudo ufw allow 22
sudo ufw allow 80
sudo ufw allow 8080
sudo ufw allow 443
sudo ufw enable


localhost/loopback

sudo iptables -t nat -I OUTPUT -p tcp -d 127.0.0.1 --dport 80 -j REDIRECT --to-ports 3000

external

sudo iptables -t nat -I PREROUTING -p tcp --dport 80 -j REDIRECT --to-ports 3000


screen -S name
screen -r name
ctr a + d   -> to detach
pkill screen





Spot : 151 Psi : 160 Phi : -92

Spot : 167 Psi : 154 Phi : -105