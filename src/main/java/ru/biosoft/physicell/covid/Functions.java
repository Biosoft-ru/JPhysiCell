package ru.biosoft.physicell.covid;

import java.util.Set;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.standard.StandardElasticContact;

public class Functions
{
    

    Cell check_for_live_neighbor_for_interaction(Cell pAttacker, double dt)
    {
        //        Set<Cell> nearby = pAttacker.cells_in_my_container();
        for( Cell cell : pAttacker.cells_in_my_container() )
        {
            if( !cell.equals( pAttacker ) && !cell.phenotype.death.dead )
                return cell;
        }
        //        int i = 0; 
        //        while( i < nearby.size() )
        //        {
        //            // don't try to kill yourself 
        //            if( nearby[i] != pAttacker && nearby[i].phenotype.death.dead == false )
        //            { return nearby[i]; }
        //            i++; 
        //        }
        return null;
    }

    Cell check_for_dead_neighbor_for_interaction(Cell pAttacker, double dt)
    {
        for( Cell cell : pAttacker.cells_in_my_container() )
        {
            if( !cell.equals( pAttacker ) && cell.phenotype.death.dead )
                return cell;
        }
        return null;
        //        Set<Cell> nearby = pAttacker.cells_in_my_container();
        //        int i = 0; 
        //        while( i < nearby.size() )
        //        {
        //            // don't try to kill yourself 
        //            if( nearby[i] != pAttacker && nearby[i].phenotype.death.dead == true )
        //            { return nearby[i]; }
        //            i++; 
        //        }
        //        return null; 
    }

    


}
