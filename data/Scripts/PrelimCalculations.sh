#!/bin/bash
if [ "$#" -ne 6 ];then
      echo "
# usage ./PrelimCalc.sh TypeOfExecution NumberOfCPU NumberOfHours/NbJobPar PathToLAMMPSbin PathToLAMMPSScripts SymmGB 
# - TypeOfExecution: seq(sequential) slurm(slurm scheduler) gnupar(gnu parallel execution)
# - NumberOfCPU: number of CPU to use for running lammps
# - PathToLAMMPSbin: path to the lammps executable
# - PathToLAMMPSScript: path to the lammps scripts (it should contain 3 files)
# - NumberOfHours/NbJobPar: either number of hours for running each slurm job (one corresponding to one free surface combination) or number of parallel job (when using gnupar option)
# - SymmGB should 0 (not a symmetric GB) or 1 (a symmetric GB)
# to take care :
# 	- well provide the path to lammps scripts
#	- well provide the good potential files in the lammps scripts
#	- take care that all atomic files are well generated in InitConfigs directory
"
	exit
fi

typeexec=$1
nbprocs=$2
nbhours=$3
path2lammpsbin=$4
path2lammpsscripts=$5
SymmGB=$6

if [ "$typeexec" = "slurm" ];then
	### MODIFY VALUES BELOW ACCORDING TO YOUR NEEDS  ###
	#SBATCH --nodes=1
	#SBATCH --ntasks-per-node=12
	#SBATCH --mem=5G
	#SBATCH --time=0:30:00
	#SBATCH --job-name=PrelimCalcs_TDOF
	
	### Load the module corresponding to the LAMMPS version you want to use
	### (type "module avail" to list available modules and versions)
	module load lammps/2Aug2023
fi

if [ -d "Bulk1" ];then 
	rm -r Bulk1
fi
if [ -d "Bulk2" ];then 
	rm -r Bulk2
fi
if [ -d "Surfaces" ];then
	rm -r Surfaces
fi
if [ -d "RelaxGB" ];then
	rm -r RelaxGB
fi

rootdir=$(pwd)

nbpup=$(ls -1 InitConfigs/UpperGrain_*.lmp | wc -l)
nbplo=$(ls -1 InitConfigs/LowerGrain_*.lmp | wc -l)
nduplo=$(awk '{print $1}' InitConfigs/NDup.dat)
ndupup=$(awk '{print $2}' InitConfigs/NDup.dat)

echo "$nbpup $nbplo" > PrelimData.dat

echo "$nbpup upper and $nbplo free surface systems found"
nbpup=$(echo $nbpup | awk '{print $1-1}')
nbplo=$(echo $nbplo | awk '{print $1-1}')
echo "Relaxing bulks.."

# First relax bulk grains
mkdir Bulk1
cp InitConfigs/LowerGrainBulk.lmp Bulk1/Bulk.lmp
cd Bulk1

if [ "$typeexec" = "slurm" ];then
	mpiexec $path2lammpsbin < $path2lammpsscripts/RelaxBulk.lin > log.lammps
else
	mpirun -n $nbprocs $path2lammpsbin -in $path2lammpsscripts/RelaxBulk.lin > log.lammps
fi

e_bulk_lo=`awk '{print $1}' EnergieBulk.dat`
n_bulk_lo=`awk '{print $2}' EnergieBulk.dat`
cd ../
mkdir Bulk2
cp InitConfigs/UpperGrainBulk.lmp Bulk2/Bulk.lmp
cd Bulk2

if [ "$typeexec" = "slurm" ];then
	mpiexec $path2lammpsbin < $path2lammpsscripts/RelaxBulk.lin > log.lammps
else
	mpirun -n $nbprocs $path2lammpsbin -in $path2lammpsscripts/RelaxBulk.lin > log.lammps
fi

e_bulk_up=`awk '{print $1}' EnergieBulk.dat`
n_bulk_up=`awk '{print $2}' EnergieBulk.dat`
cd ../

echo "${e_bulk_up} ${n_bulk_up} ${e_bulk_lo} ${n_bulk_lo}" >> PrelimData.dat

echo "done"
echo "Computing upper surface energies.."

# Compute surfaces energies
mkdir Surfaces

mkdir Surfaces/Upper
awk -v eb1=${e_bulk_up} -v nb1=${n_bulk_up} '{if($2=="MyVars_eg") print $1,$2,$3,eb1;else if($2=="MyVars_ng") print $1,$2,$3,nb1;else print $0}' ${path2lammpsscripts}/RelaxSurface.lin > Surfaces/Upper/RelaxSurface.lin

mkdir Surfaces/Lower
awk -v eb1=${e_bulk_lo} -v nb1=${n_bulk_lo} '{if($2=="MyVars_eg") print $1,$2,$3,eb1;else if($2=="MyVars_ng") print $1,$2,$3,nb1;else print $0}' ${path2lammpsscripts}/RelaxSurface.lin > Surfaces/Lower/RelaxSurface.lin

declare -a gamma_up
declare -a nbg_up
cd Surfaces/Upper

for p in $(seq 0 1 $nbpup);do
	mkdir Surf_${p}
	cd Surf_${p}
	cp ${rootdir}/InitConfigs/UpperGrain_${p}.lmp Surface.lmp
	nbg_up[$p]=`awk -v n=$ndupup 'NR==3{print $1*n}' Surface.lmp` # TODO warning here if the format of file changes
	if [ "$typeexec" = "slurm" ];then
		mpiexec $path2lammpsbin < ../RelaxSurface.lin > log.lammps
	elif [ "$typeexec" = "seq" ];then
		mpirun -n $nbprocs $path2lammpsbin -in ../RelaxSurface.lin > log.lammps
	fi
	if [ "$typeexec" != "gnupar" ];then
		gamma_up[$p]=`awk '{print $1}' gamma.dat`
		echo "$p ${nbg_up[$p]} ${gamma_up[$p]}" >> ${rootdir}/PrelimData.dat
		file=$( ls -v Relaxed_*.cfg | tail -1)
		mv $file tmp
		rm *.cfg
		rm Surface.lmp
		mv tmp RelaxedSurface.cfg
	fi
	cd ../
done
if [ "$typeexec" = "gnupar" ];then
	export nbprocs
	export path2lammpsbin
	export rootdir
	parallel -j ${nbhours} --env nbprocs,path2lammpsbin,rootdir '
		p={1}	
		cd Surf_$p
		nbg=$(awk '\''NR==3{print $1}'\'' Surface.lmp) # TODO warning here if the format of file changes
		mpirun -n $nbprocs $path2lammpsbin -in ../RelaxSurface.lin > log.lammps
		gamma=$(awk '\''{print $1}'\'' gamma.dat)
		echo "up $p $nbg $gamma" >> ${rootdir}/PrelimData.dat
		file=$( ls -v Relaxed_*.cfg | tail -1)
		mv $file tmp
		rm *.cfg
		rm Surface.lmp
		mv tmp RelaxedSurface.cfg
		cd ../
	' ::: $(seq 0 1 $nbpup)
	for i in $(seq 0 1 $nbpup);do
		nbg_up[$i]=$(awk -v i=$i -v n=$ndupup '{if(($1=="up")&&($2==i)) print $3*n}' ${rootdir}/PrelimData.dat )
		gamma_up[$i]=$(awk -v i=$i '{if(($1=="up")&&($2==i)) print $4}' ${rootdir}/PrelimData.dat )
	done
fi
cd ../

echo "done"
echo "Computing lower surface energies.."

declare -a gamma_lo
declare -a nbg_lo
cd Lower

for p in $(seq 0 1 $nbplo);do
	mkdir Surf_${p}
	cd Surf_${p}
	cp ${rootdir}/InitConfigs/LowerGrain_${p}.lmp Surface.lmp
	nbg_lo[$p]=`awk -v n=$nduplo 'NR==3{print $1*n}' Surface.lmp` # TODO warning here if the format of file changes
	if [ "$typeexec" = "slurm" ];then
		mpiexec $path2lammpsbin < ../RelaxSurface.lin > log.lammps
	elif [ "$typeexec" = "seq" ];then
		mpirun -n $nbprocs $path2lammpsbin -in ../RelaxSurface.lin > log.lammps
	fi
	if [ "$typeexec" != "gnupar" ];then
		gamma_lo[$p]=`awk '{print $1}' gamma.dat`
		echo "$p ${nbg_lo[$p]} ${gamma_lo[$p]}" >> ${rootdir}/PrelimData.dat
		file=$( ls -v Relaxed_*.cfg | tail -1)
		mv $file tmp
		rm *.cfg
		rm Surface.lmp
		mv tmp RelaxedSurface.cfg
	fi
	cd ../
done
if [ "$typeexec" = "gnupar" ];then
	export nbprocs
	export path2lammpsbin
	export rootdir
	parallel -j ${nbhours} --env nbprocs,path2lammpsbin,rootdir '
		p={1}	
		cd Surf_$p
		nbg=$(awk '\''NR==3{print $1}'\'' Surface.lmp) # TODO warning here if the format of file changes
		mpirun -n $nbprocs $path2lammpsbin -in ../RelaxSurface.lin > log.lammps
		gamma=$(awk '\''{print $1}'\'' gamma.dat)
		echo "lo $p $nbg $gamma" >> ${rootdir}/PrelimData.dat
		file=$( ls -v Relaxed_*.cfg | tail -1)
		mv $file tmp
		rm *.cfg
		rm Surface.lmp
		mv tmp RelaxedSurface.cfg
		cd ../
	' ::: $(seq 0 1 $nbplo)
	for i in $(seq 0 1 $nbplo);do
		nbg_lo[$i]=$(awk -v i=$i -v n=$nduplo '{if(($1=="lo")&&($2==i)) print $3*n}' ${rootdir}/PrelimData.dat )
		gamma_lo[$i]=$(awk -v i=$i '{if(($1=="lo")&&($2==i)) print $4}' ${rootdir}/PrelimData.dat )
	done
fi
echo "done"

cd ../../

# now generate the script for launching calculations of GB relaxation
echo "UpSurfaceIndex LoSurfaceIndex XShift YShift Gamma(J/m2)" > TotEnergies.dat
mkdir RelaxGB
cd RelaxGB
for plo in $(seq 0 1 $nbplo);do
	pup_beg=0
	if [ $SymmGB -eq 1 ];then
		pup_beg=$plo
	fi
	for pup in $(seq $pup_beg 1 $nbpup);do
		echo "Computing GB energies for upper free surface $pup and lower free surface $plo"
		mkdir GB_Low_${plo}_Up_${pup}
		cd GB_Low_${plo}_Up_${pup}
		awk -v ng1=${nbg_up[$pup]} -v ng2=${nbg_lo[$plo]} -v gam1=${gamma_up[$pup]} -v gam2=${gamma_lo[$plo]} -v eb1=$e_bulk_up -v nb1=$n_bulk_up -v eb2=$e_bulk_lo -v nb2=$n_bulk_lo '{
		if($2=="MyVars_eg1"){
			print $1,$2,$3,eb1;
		}else if($2=="MyVars_eg2"){
			print $1,$2,$3,eb2;
		}else if($2=="MyVars_nb1"){
			print $1,$2,$3,nb1;
		}else if($2=="MyVars_nb2"){
			print $1,$2,$3,nb2;
		}else if($2=="MyVars_ng1"){
			print $1,$2,$3,ng1;
		}else if($2=="MyVars_ng2"){
			print $1,$2,$3,ng2;
		}else if($2=="MyVars_gamma1"){
			print $1,$2,$3,gam1;
		}else if($2=="MyVars_gamma2"){
			print $1,$2,$3,gam2;
		}else{
		       	print $0
		}}' ${path2lammpsscripts}/RelaxGB.lin > RelaxGB.lin
		# search number of shift
		tempstr=$(ls -v ${rootdir}/InitConfigs/GB_PlaneLow_${plo}_PlaneUp_${pup}_Shift_*.lmp | tail -1)
		tempstr=${tempstr#"$rootdir"}
		nbsx=$(echo $tempstr | awk -F "_" '{print $7}')
		nbsy=$(echo $tempstr | awk -F "_" '{print $8}')
		nbsy=$(echo $nbsy | awk -F "." '{print $1}')

		if [ "$typeexec" = "slurm" ];then

		cat << EOF > LaunchCalc_${plo}_${pup}.slurm
#!/bin/bash
  
### MODIFY VALUES BELOW ACCORDING TO YOUR NEEDS  ###
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=${nbprocs}
#SBATCH --mem=5G
#SBATCH --time=${nbhours}:00:00
#SBATCH --job-name=tdof_${plo}_${pup}

### Load the module corresponding to the LAMMPS version you want to use
### (type "module avail" to list available modules and versions)
module load lammps/2Aug2023
module load gcc/11.3.1/compilers
module load cmake/3.29.2

echo "ShiftX ShiftY GammaGB(J/m2)" > Energies.dat
for sx in \$(seq 0 1 $nbsx);do
	for sy in \$(seq 0 1 $nbsy);do
		mkdir Shift_\${sx}_\${sy}
		cd Shift_\${sx}_\${sy}
		cp ${rootdir}/InitConfigs/GB_PlaneLow_${plo}_PlaneUp_${pup}_Shift_\${sx}_\${sy}.lmp GB.lmp
		mpiexec lmp < ../RelaxGB.lin > log.lammps
		gammaGB=\`awk '{print \$1}' EnergieGB.dat\`
		echo "${pup} ${plo} \${sx} \${sy} \${gammaGB}" >> ${rootdir}/TotEnergies.dat
		echo "\${sx} \${sy} \${gammaGB}" >> ../Energies.dat
		file=\$( ls -v Relaxed_*.cfg | tail -1)
		mv \$file tmp
		rm *.cfg
		rm GB.lmp
		mv tmp RelaxedGB.cfg
		cd ../
	done
done

EOF
		else
			for sx in $(seq 0 1 $nbsx);do
				for sy in $(seq 0 1 $nbsy);do
					mkdir Shift_${sx}_${sy}
					cd Shift_${sx}_${sy}
					cp ${rootdir}/InitConfigs/GB_PlaneLow_${plo}_PlaneUp_${pup}_Shift_${sx}_${sy}.lmp GB.lmp
					if [ "$typeexec" = "seq" ];then
						mpirun -n $nbprocs $path2lammpsbin -in ../RelaxGB.lin > log.lammps
						gammaGB=`awk '{print $1}' EnergieGB.dat`
						echo "${pup} ${plo} ${sx} ${sy} ${gammaGB}" >> ${rootdir}/TotEnergies.dat
						echo "${sx} ${sy} ${gammaGB}" >> ../Energies.dat
						file=$( ls -v Relaxed_*.cfg | tail -1)
						mv $file tmp
						rm *.cfg
						rm GB.lmp
						mv tmp RelaxedGB.cfg
					fi
					cd ../
				done
			done
		fi
		if [ "$typeexec" = "gnupar" ];then
			export nbprocs
			export path2lammpsbin
			export rootdir
			export pup
			export plo
			parallel -j ${nbhours} --env nbprocs,path2lammpsbin,rootdir,pup,plo '
				sx={1}	
				sy={2}	
				cd Shift_${sx}_${sy}
				mpirun -n $nbprocs $path2lammpsbin -in ../RelaxGB.lin > log.lammps
				gammaGB=$(awk '\''{print $1}'\'' EnergieGB.dat)
				echo "${pup} ${plo} ${sx} ${sy} ${gammaGB}" >> ${rootdir}/TotEnergies.dat
				echo "${sx} ${sy} ${gammaGB}" >> ../Energies.dat
				file=$( ls -v Relaxed_*.cfg | tail -1)
				mv $file tmp
				rm *.cfg
				rm GB.lmp
				mv tmp RelaxedGB.cfg
				cd ../
			' ::: $(seq 0 1 $nbsx) ::: $(seq 0 1 $nbsy)
		fi
		cd ../
		echo "done"
	done
done

cd ../

# generate last shell script for launching sims
if [ "$typeexec" = "slurm" ];then
	cat << EOF > SubmitJobs_tdof.sh
cd RelaxGB
for i in \`ls -d GB_*/\`;do
	cd \$i
	sbatch LaunchCalc_*.slurm
	cd ../
done
EOF
fi
