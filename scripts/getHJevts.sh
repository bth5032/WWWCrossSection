#!/bin/bash

function makeERL {
  makeERL_out=`ssh uaf-10 "cat /home/users/haweber/Tribosons/logs/$1" | awk -F ':' '{print "make_tuple("$3","$1","$2"),"}' | xargs`
  echo -n ${makeERL_out::${#makeERL_out}-1}
}

if [[ $1 == "all" ]]
then
  echo "//---------"
  echo "// EE"
  echo "//---------"
  echo "//tt1l"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListEE_tt1l.log`"};"
  echo "//WZ"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListEE_WZ.log`"};"
  echo "//Wjets"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListEE_Wjets.log`"};"
  echo "//WWW"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListEE_WWW.log`"};"
  echo ""
  echo ""
  echo "//---------"
  echo "// MM"
  echo "//---------"
  echo "//tt1l"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListMM_tt1l.log`"};"
  echo "//WZ"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListMM_WZ.log`"};"
  echo "//Wjets"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListMM_Wjets.log`"};"
  echo "//WWW"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListMM_WWW.log`"};"
  echo "//---------"
  echo "// EM"
  echo "//---------"
  echo "//tt1l"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListEM_tt1l.log`"};"
  echo "//WZ"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListEM_WZ.log`"};"
  echo "//Wjets"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListEM_Wjets.log`"};"
  echo "//WWW"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventListEM_WWW.log`"};"
  echo "//---------"
  echo "// 0SFOS"
  echo "//---------"
  echo "//tt1l"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList0SFOS_tt1l.log`"};"
  echo "//WZ"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList0SFOS_WZ.log`"};"
  echo "//Wjets"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList0SFOS_Wjets.log`"};"
  echo "//WWW"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList0SFOS_WWW.log`"};"
  echo "//---------"
  echo "// 1SFOS"
  echo "//---------"
  echo "//tt1l"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList1SFOS_tt1l.log`"};"
  echo "//WZ"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList1SFOS_WZ.log`"};"
  echo "//Wjets"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList1SFOS_Wjets.log`"};"
  echo "//WWW"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList1SFOS_WWW.log`"};"
  echo "//---------"
  echo "// 2SFOS"
  echo "//---------"
  echo "//tt1l"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList2SFOS_tt1l.log`"};"
  echo "//WZ"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList2SFOS_WZ.log`"};"
  echo "//Wjets"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList2SFOS_Wjets.log`"};"
  echo "//WWW"
  echo "//set<tuple<long,long,long>>  inspection_set_erl = {"`makeERL eventList2SFOS_WWW.log`"};"
elif [[ $# < 2 ]]
then
  echo "Usage: getHJevts <SR> <Sample>  OR getHJevts all"
  echo "-------------"
  echo "SRs: "
  echo "-------------"
  echo "0SFOS"
  echo "1SFOS"
  echo "2SFOS"
  echo "EE"
  echo "MM"
  echo "EM"
  echo "-------------"
  echo "Recommended Samples: "
  echo "-------------"
  echo "tt1l"
  echo "WZ"
  echo "WWW"
  echo "WJets"
  echo "-------------"
  echo "Everything in Dir: "
  echo "-------------"
  ssh uaf-10 "ls /home/users/haweber/Tribosons/logs/eventList0SFOS*" | sed 's/^.*eventList0SFOS_\(.*\)\.log$/\1/g'
  exit 0
else
  makeERL eventList$1_$2.log
fi