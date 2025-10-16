## Burnup schemes : 1L or 2L flux calculations
# Date : 25/09/2025

#
import lifo
import cle2000

def burnup2LScheme(lib_lcm,TRACK,TF_EXC,TRACK_L1,TF_EXC_L1,TRACK_SSH,TF_EXC_SSH,StepList,name_compo,ssh_solver,lvl1_solver,name_geom):
    """
    PARAMETER COMPO LIBRARY TRACK TF_EXC 
    TRACK_1L TF_EXC_L1 TRACK_SS TF_EXC_SS StepList ::
    ::: LINKED_LIST COMPO ;
    ::: LINKED_LIST LIBRARY ; 
    ::: LINKED_LIST TRACK ;
    ::: SEQ_BINARY TF_EXC ;
    ::: LINKED_LIST TRACK_1L ;
    ::: SEQ_BINARY TF_EXC_L1 ;
    ::: LINKED_LIST TRACK_SS ;
    ::: SEQ_BINARY TF_EXC_SS ;
    ::: LINKED_LIST StepList ; ;

    STRING name_compo ssh_sol lvl1_sol ;
    :: >>name_compo<< >>ssh_sol<< >>lvl1_sol<< ;

    """
    print(f"inputs : ssh_solver = {ssh_solver}, lvl1_solver = {lvl1_solver}")

    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty("COMPO","LCM")
    myLifo.push(lib_lcm)
    myLifo.push(TRACK)
    myLifo.push(TF_EXC)
    myLifo.push(TRACK_L1)
    myLifo.push(TF_EXC_L1)
    myLifo.push(TRACK_SSH)
    myLifo.push(TF_EXC_SSH)
    myLifo.push(StepList)
    myLifo.push(name_compo)
    myLifo.push(ssh_solver)
    myLifo.push(lvl1_solver)
    myLifo.push(name_geom)
    myLifo.lib()

    # Execution 
    BWR_burnup_2L = cle2000.new('BU_A_2L',myLifo,1)
    BWR_burnup_2L.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node("COMPO")
    return pyCOMPO

def burnup1LScheme(lib_lcm,TRACK,TF_EXC,TRACK_SSH,TF_EXC_SSH,StepList,name_compo,ssh_solver,name_geom):
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
    myLifo.lib()

    # Execution 
    BWR_burnup_1L = cle2000.new('BU_A_1L',myLifo,1)
    BWR_burnup_1L.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node("COMPO")
    return pyCOMPO
