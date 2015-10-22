echo "You'll need Python 2, nose and"
echo "MDAnalysis (and a few other" 
echo "dependencies!). Please run from"
echo "the ligand_interaction_network"
echo "root directory."
echo "______________________________"
echo "Here we go..."
echo ""
rm -Rf deleteme_*
python lintools/extract_data.py -j deleteme_shelltest -i data/docked_ami_100ns_nowater.gro -t data/docked_ami_50ns_1_100_whole_no_water_skip100.xtc
python lintools
python lintools/visualisation.py -j deleteme_shelltest
nosetests tests/tests.py
rm -Rf deleteme_*
echo ""
echo "That should have worked..."