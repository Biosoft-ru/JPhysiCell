package ru.biosoft.physicell.covid;

import java.util.Set;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;

public class CD8_Tcell_mechanics extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt) throws Exception
    {
        //    void CD8_Tcell_mechanics( Cell pCell, Phenotype phenotype, double dt )
        //    {
        int debris_index = pCell.getMicroenvironment().findDensityIndex( "debris" );

        if( phenotype.death.dead == true )
        {
            pCell.functions.updatePhenotype = null;
            pCell.functions.customCellRule = null;

            phenotype.secretion.secretionRates[debris_index] = pCell.customData.get( "debris_secretion_rate" );
            return;
        }

        // bounds check 
        if( ModelCovid.check_for_out_of_bounds( pCell, 10.0 ) )
        {
            //            #pragma omp critical
            //                    {
            ( (ModelCovid)pCell.getModel() ).cells_to_move_from_edge.add( pCell );
            //                    }
            // replace_out_of_bounds_cell( pCell, 10.0 );
            // return; 
        }

        // if I am not adhered to a cell, turn motility on 
        if( pCell.state.neighbors.size() == 0 )
        {
            phenotype.motility.isMotile = true;
        }
        else
        {
            phenotype.motility.isMotile = false;
        }

        // check for contact with infected cell 

        // if I'm adhered to something ... 
        if( pCell.state.numberAttachedCells() > 0 ) // pCell.state.neighbors.size() > 0 )
        {
            // decide whether to detach 
            boolean detach_me = false;

            if( pCell.getModel().getRNG().UniformRandom() < dt / ( pCell.customData.get( "cell_attachment_lifetime" ) + 1e-15 ) )
            {
                detach_me = true;
            }

            // if I detach, go through the process 
            if( detach_me )
            {
                pCell.removeAllAttachedCells();//.remove_all_attached_cells(); 
                // resume motile behavior 
                phenotype.motility.isMotile = true;
            }
            return;
        }

        // I'm not attached, look for cells nearby and try to attach

        // if this returns non-null, we're now attached to a cell 
        if( immune_cell_check_neighbors_for_attachment( pCell, dt ) != null )
        {
            // set motility off 
            phenotype.motility.isMotile = false;
            return;
        }
        phenotype.motility.isMotile = true; // I suggest eliminating this. 
    }
    
    static boolean attempt_immune_cell_attachment(Cell pAttacker, Cell pTarget, double dt)
    {
        // if the target is not infected, give up 
        if( pTarget.customData.get( "viral_protein" ) < pAttacker.customData.get( "TCell_detection" ) )
        {
            return false;
        }

        // if the target cell is dead, give up 
        if( pTarget.phenotype.death.dead == true )
        {
            return false;
        }

        // if the target cell is too far away, give up 
        double[] displacement = VectorUtil.newDiff( pTarget.position, pAttacker.position );
        double distance_scale = VectorUtil.norm( displacement );

        // better: use mechanics constants 
        if( distance_scale > pAttacker.customData.get( "max_attachment_distance" ) )
        {
            return false;
        }

        // now, get the attachment probability 

        double attachment_probability = pAttacker.customData.get( "cell_attachment_rate" ) * dt;

        // don't need to cap it at 1.00: if prob > 100%, 
        // then this statement always evaluates as true, 
        // just the same as capping probability at 100% 
        if( pAttacker.getModel().getRNG().UniformRandom() <= attachment_probability )
        {
            Cell.attachcCells( pAttacker, pTarget );
            return true;
        }

        return false;
    }

    static Cell immune_cell_check_neighbors_for_attachment(Cell pAttacker, double dt)
    {
        Set<Cell> nearby = pAttacker.cells_in_my_container();
        for( Cell cell : nearby )
        {
            if( !cell.equals( pAttacker ) )
            {
                if( attempt_immune_cell_attachment( pAttacker, cell, dt ) )
                {
                    return cell;
                }
            }
            //              i++; 
        }
        //        int i = 0; 
        //        while( i < nearby.size() )
        //        {
        //            // don't try to kill yourself 
        //            if( nearby[i] != pAttacker )
        //            {
        //                if( attempt_immune_cell_attachment( pAttacker, nearby[i] , dt ) )
        //                {
        //                    return nearby[i]; 
        //                }
        //            }
        //            i++; 
        //        }

        return null;
    }
}