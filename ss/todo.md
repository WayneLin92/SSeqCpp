# Main
* add Z to ss
* add tmf_C2_Ceta_Cnu to ss
* add CW_2_V_eta = RP3_5
* Adams scheduler support d2
* plot extensions
* check bottom cell status
 

# Adams
* Brown-Gitler
* Fast multiplication using extensions

# ss
* Add more cofseq
* Improve rename_gen
* smaller log.db
* Cache contradictions

* Cache multiplications
* Simplify pi generators and relations
* Fast sync
* Test differential monomial orderings for mod
* Optimize deduce diff for dr=0 and increase r.
* Support stem_trunc

# html
* multiplication by theta and g
* Mark differentials that cannot be zero

# Others
## tasks
./Adams export_mod C2nu S0 200
./Adams export_mod Ctheta5sq S0 200
./Adams export_mod C2sigma S0 200
./Adams export_mod Cetasigma S0 200
./Adams export_map C2nu CP1_4 200
./Adams export_map CW_2_eta_theta5sq Ctheta5sq 200
./Adams export_map CW_nu_eta_2_etasigma Cetasigma 200
./Adams export_map CW_sigma_2sigma C2sigma 200
./Adams export_map C2sigma CW_2sigma_sigma 200
./Adams export_map Cetasigma CW_etasigma_2_eta_nu 200
./Adams export_map Ctheta5sq CW_theta5sq_eta_2 200
./Adams export_map Ctheta5sq C2_by_eta 200


cp C2nu_AdamsSS.db AdamsSS
cp Ctheta5sq_AdamsSS.db AdamsSS
cp C2sigma_AdamsSS.db AdamsSS
cp Cetasigma_AdamsSS.db AdamsSS

cp map_AdamsSS_C2nu__CP1_4.db AdamsSS
cp map_AdamsSS_CW_2_eta_theta5sq__Ctheta5sq.db AdamsSS
cp map_AdamsSS_CW_nu_eta_2_etasigma__Cetasigma.db AdamsSS
cp map_AdamsSS_CW_sigma_2sigma__C2sigma.db AdamsSS
cp map_AdamsSS_C2sigma__CW_2sigma_sigma.db AdamsSS
cp map_AdamsSS_Cetasigma__CW_etasigma_2_eta_nu.db AdamsSS
cp map_AdamsSS_Ctheta5sq__CW_theta5sq_eta_2.db AdamsSS
cp map_AdamsSS_Ctheta5sq__C2_by_eta.db AdamsSS

./Adams export_mod C4 S0 200
./Adams export_mod C8 S0 200
./Adams export_mod Cetasq S0 200
./Adams export_mod Cetacube S0 200
./Adams export_map C4 S0 200
./Adams export_map C8 S0 200
./Adams export_map Cetasq S0 200
./Adams export_map Cetacube S0 200

cp C4_AdamsSS.db AdamsSS
cp C8_AdamsSS.db AdamsSS
cp Cetasq_AdamsSS.db AdamsSS
cp Cetacube_AdamsSS.db AdamsSS
cp map_AdamsSS_C4__S0.db AdamsSS
cp map_AdamsSS_C8__S0.db AdamsSS
cp map_AdamsSS_Cetasq__S0.db AdamsSS
cp map_AdamsSS_Cetacube__S0.db AdamsSS