conda activate moose
conda update --all
cd ~/projects/moose
git fetch origin
git reset --hard
git rebase origin/master
cd ~/projects/moose_beh
make clobberall
make -j16
./run_tests -j16
