package ru.biosoft.physicell.sample_projects.cancer_immune;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.custom_cell_rule;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellUtilities;

public class ImmuneCellRule implements custom_cell_rule
{
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        int attach_lifetime_i = pCell.custom_data.findVariableIndex( "attachment_lifetime" );

        if( phenotype.death.dead == true )
        {
            // the cell death functions don't automatically turn off custom functions, 
            // since those are part of mechanics. 

            // Let's just fully disable now. 
            pCell.functions.custom_cell_rule = null;
            return;
        }

        // if I'm docked
        if( pCell.state.numberAttachedCells() > 0 )
        {
            // attempt to kill my attached cell
            Cell attached = pCell.state.attachedCells.iterator().next();///[0];
            boolean detach_me = false;

            if( immune_cell_attempt_apoptosis( pCell, attached, dt ) )
            {
                immune_cell_trigger_apoptosis( pCell, attached );
                detach_me = true;
            }

            // decide whether to detach 
            if( PhysiCellUtilities.UniformRandom() < dt / ( pCell.custom_data.get( attach_lifetime_i ) + 1e-15 ) )
            {
                detach_me = true;
            }

            // if I dettach, resume motile behavior 
            if( detach_me )
            {
                Cell.detachCells( pCell, attached );
                phenotype.motility.isMotile = true;
            }
            return;
        }

        // I'm not docked, look for cells nearby and try to docked
        // if this returns non-NULL, we're now attached to a cell 
        if( immune_cell_check_neighbors_for_attachment( pCell, dt ) != null )
        {
            // set motility off 
            phenotype.motility.isMotile = false;
            return;
        }
        phenotype.motility.isMotile = true;
    }

    static boolean immune_cell_trigger_apoptosis(Cell pAttacker, Cell pTarget)
    {
        int apoptosis_model_index = pTarget.phenotype.death.findDeathModelIndex( "apoptosis" );
        if( pTarget.phenotype.death.dead == true )
        {
            return false; // if the Target cell is already dead, don't bother!
        }
        pTarget.startDeath( apoptosis_model_index );
        return true;
    }

    static boolean immune_cell_attempt_apoptosis(Cell pAttacker, Cell pTarget, double dt)
    {
        int oncoprotein_i = pTarget.custom_data.findVariableIndex( "oncoprotein" );
        int apoptosis_model_index = pTarget.phenotype.death.findDeathModelIndex( "apoptosis" );
        int kill_rate_index = pAttacker.custom_data.findVariableIndex( "kill_rate" );

        double oncoprotein_saturation = pAttacker.custom_data.get( "oncoprotein_saturation" ); // 2.0; 
        double oncoprotein_threshold = pAttacker.custom_data.get( "oncoprotein_threshold" ); // 0.5; // 0.1; 
        double oncoprotein_difference = oncoprotein_saturation - oncoprotein_threshold;

        if( pTarget.custom_data.get( oncoprotein_i ) < oncoprotein_threshold )
        {
            return false;
        }

        double scale = pTarget.custom_data.get( oncoprotein_i );
        scale -= oncoprotein_threshold;
        scale /= oncoprotein_difference;
        if( scale > 1.0 )
        {
            scale = 1.0;
        }

        if( PhysiCellUtilities.UniformRandom() < pAttacker.custom_data.get( kill_rate_index ) * scale * dt )
        {
            //          std::cout << "\t\t kill!" << " " << pTarget.custom_data[oncoprotein_i] << std::endl; 
            return true;
        }
        return false;
    }

    public static Cell immune_cell_check_neighbors_for_attachment(Cell pAttacker, double dt)
    {
        for( Cell nearbyCell : pAttacker.cells_in_my_container() )
        {
            if( nearbyCell != pAttacker )// don't try to kill yourself 
            {
                if( immune_cell_attempt_attachment( pAttacker, nearbyCell, dt ) )
                {
                    return nearbyCell;
                }
            }
        }
        return null;
    }

    static boolean immune_cell_attempt_attachment(Cell pAttacker, Cell pTarget, double dt)
    {
        int oncoprotein_i = pTarget.custom_data.findVariableIndex( "oncoprotein" );
        int attach_rate_i = pAttacker.custom_data.findVariableIndex( "attachment_rate" );

        double oncoprotein_saturation = pAttacker.custom_data.get( "oncoprotein_saturation" );
        double oncoprotein_threshold = pAttacker.custom_data.get( "oncoprotein_threshold" );
        double oncoprotein_difference = oncoprotein_saturation - oncoprotein_threshold;

        double max_attachment_distance = pAttacker.custom_data.get( "max_attachment_distance" );
        double min_attachment_distance = pAttacker.custom_data.get( "min_attachment_distance" );
        double attachment_difference = max_attachment_distance - min_attachment_distance;

        if( pTarget.custom_data.get( oncoprotein_i ) > oncoprotein_threshold && pTarget.phenotype.death.dead == false )
        {
            double[] displacement = VectorUtil.newDiff( pTarget.position, pAttacker.position );
            double distance_scale = VectorUtil.norm( displacement );
            if( distance_scale > max_attachment_distance )
            {
                return false;
            }

            double scale = pTarget.custom_data.get( oncoprotein_i );
            scale -= oncoprotein_threshold;
            scale /= oncoprotein_difference;
            if( scale > 1.0 )
            {
                scale = 1.0;
            }

            distance_scale *= -1.0;
            distance_scale += max_attachment_distance;
            distance_scale /= attachment_difference;
            if( distance_scale > 1.0 )
            {
                distance_scale = 1.0;
            }

            if( PhysiCellUtilities.UniformRandom() < pAttacker.custom_data.get( attach_rate_i ) * scale * dt * distance_scale )
            {
                //              std::cout << "\t attach!" << " " << pTarget.custom_data[oncoprotein_i] << std::endl; 
                Cell.attachcCells( pAttacker, pTarget );
            }
            return true;
        }
        return false;
    }
}