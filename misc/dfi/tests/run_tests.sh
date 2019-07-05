module load psi4/1.2.1.py36_oepdev-1.0.2 
psi4 --psiapi-path
python -c "import psi4"

python DFI_JK.py
python DFI_J.py
