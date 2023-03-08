# Re-Glyco

The Python software utilizes molecular dynamics (MD) simulation results from the glycoshape database to re-glycosylate protein structures predicted by alphafold, a deep learning method for protein structure prediction. The input to the software is a protein structure output by alphafold, and the output is a modified protein structure with glycans (sugar molecules) added at appropriate sites. The incorporation of glycans is achieved through the use of MD simulation results, which ensures that the resulting glycosylated protein structure is physically realistic and stable. This software can be useful for researchers studying the role of glycans in protein function and the impact of glycosylation on protein structure and stability.


# Installation
```
conda create -n reglyco python=3.10
conda activate reglyco
pip install -r requirements.txt
streamlit run main.py
```


To implement!
1. convegence thing with max iter limit.
2. warning user that we used non standard angle to attach, as standard failed.
3. Implement Cluster.

python3 -m pip -r requirements.txt
export PATH="$HOME/.local/bin:$PATH"


sudo iptables -P INPUT ACCEPT
sudo iptables -P OUTPUT ACCEPT
sudo iptables -P FORWARD ACCEPT
sudo iptables -F
sudo ufw 22
sudo ufw 80
sudo ufw 8080
sudo ufw 443
sudo ufw enable


localhost/loopback

sudo iptables -t nat -I OUTPUT -p tcp -d 127.0.0.1 --dport 80 -j REDIRECT --to-ports 3000

external

sudo iptables -t nat -I PREROUTING -p tcp --dport 80 -j REDIRECT --to-ports 3000

screen -S name
screen -r name
ctr a + d   -> to detach
pkill screen
