
function getTables {
  REGION=$1
  ssh uaf-10 "source ~/.bash_profile> /dev/null 2>&1  && rootenv > /dev/null 2>&1  && ev > /dev/null 2>&1 && cd Projects/GIT/WWWCrossSection > /dev/null 2>&1  && echo '\subsection*{All Bgs}' && python scripts/makeFRTables.py --bg -s $REGION && echo && echo && echo '\subsection*{\"Fake\" Backgrounds}' && python scripts/makeFRTables.py --frbg -s $REGION && echo && echo && echo '\subsection*{Just WWW}' && python scripts/makeFRTables.py --www -s $REGION && echo && echo && echo '\subsection*{Everything}' && python scripts/makeFRTables.py --all -s $REGION && echo && echo '\\newpage' && echo "| pbcopy
}

if [[ $1 == "all" ]]
then
  echo "Doing Baseline: "
  getTables Baseline
  read -p "Press enter to continue"

  echo "Doing Mai: "
  getTables Baseline
  read -p "Press enter to continue"

  echo "Doing TightIso: "
  getTables TightIso
  read -p "Press enter to continue"

  echo "Doing LooseIso: "
  getTables LooseIso
  read -p "Press enter to continue"

  echo "Doing Btags: "
  getTables Btags
  read -p "Press enter to continue"

  echo "Doing JetMass: "
  getTables JetMass

else  
  echo "Running getTables "$1
  getTables $1
fi