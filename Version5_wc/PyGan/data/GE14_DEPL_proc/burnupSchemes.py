## Burnup schemes : 1L or 2L flux calculations
# Date : 25/09/2025

#
import lifo
import cle2000

def burnup2LScheme(lib_lcm,trk_ssh,trkfil_ssh,trk_l1,trkfil_l1,trk_l2,trkfil_l2,StepList,name_compo,assembly_id,specific_power,rates_normalisation_option,max_sph_gr):
    """
    * ---
    * INPUT & OUTPUT PARAMETERS
    * ---
    PARAMETER COMPO LIBRARY TRKSSH TRKFILSSH 
        TRKFLUXL1 TRKFILFLUXL1 TRKFLUXL2 TRKFILFLUXL2 StepList ::
    ::: LINKED_LIST COMPO ;
    ::: LINKED_LIST LIBRARY ; 
    ::: LINKED_LIST TRKSSH ;
    ::: SEQ_BINARY TRKFILSSH ;
    ::: LINKED_LIST TRKFLUXL1 ;
    ::: SEQ_BINARY TRKFILFLUXL1 ;
    ::: LINKED_LIST TRKFLUXL2 ;
    ::: SEQ_BINARY TRKFILFLUXL2 ;
    ::: LINKED_LIST StepList ; ;

    STRING name_compo assbly_id ;
    :: >>name_compo< >>assbly_id<< ;
    DOUBLE Dspec_pow ;
    :: >>Dspec_pow<< ; ! SPECIFIC POWER MW/t
    STRING rates_norm ;
    :: >>rates_norm<< ; ! OPTION FOR RATES NORMALISATION: "NONE", "QFIS" (OpenMC), "EDP0" (Serpent)

    """

    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty("COMPO","LCM")
    myLifo.push(lib_lcm)
    myLifo.push(trk_ssh)
    myLifo.push(trkfil_ssh)
    myLifo.push(trk_l1)
    myLifo.push(trkfil_l1)
    myLifo.push(trk_l2)
    myLifo.push(trkfil_l2)
    myLifo.push(StepList)
    myLifo.push(name_compo)
    myLifo.push(assembly_id)
    myLifo.push(specific_power)
    myLifo.push(rates_normalisation_option)
    myLifo.push(max_sph_gr)
    myLifo.lib()

    # Execution 
    BWR_burnup_2L = cle2000.new('BU_A_2L',myLifo,1)
    BWR_burnup_2L.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node("COMPO")
    return pyCOMPO

def burnup1LScheme(lib_lcm,TRACK,TF_EXC,TRACK_SSH,TF_EXC_SSH,StepList,name_compo,ssh_solver,name_geom, specific_power, rates_normalisation_option):
    """
    PARAMETER COMPO LIBRARY TRACK TF_EXC TRACK_SS TF_EXC_SS StepList ::
    ::: LINKED_LIST COMPO ;
    ::: LINKED_LIST LIBRARY ; 
    ::: LINKED_LIST TRACK ;
    ::: SEQ_BINARY TF_EXC ;
    ::: LINKED_LIST TRACK_SS ;
    ::: SEQ_BINARY TF_EXC_SS ;
    ::: LINKED_LIST StepList ; ;

    STRING name_compo ssh_sol ;
    :: >>name_compo<< >>ssh_sol<< ;

    """
    print(f"inputs : ssh_solver = {ssh_solver}")

    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty("COMPO","LCM")
    myLifo.push(lib_lcm)
    myLifo.push(TRACK)
    myLifo.push(TF_EXC)
    myLifo.push(TRACK_SSH)
    myLifo.push(TF_EXC_SSH)
    myLifo.push(StepList)
    myLifo.push(name_compo)
    myLifo.push(ssh_solver)
    myLifo.push(name_geom)
    myLifo.push(specific_power)
    myLifo.push(rates_normalisation_option)
    myLifo.lib()

    # Execution 
    BWR_burnup_1L = cle2000.new('BU_A_1L',myLifo,1)
    BWR_burnup_1L.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node("COMPO")
    return pyCOMPO
