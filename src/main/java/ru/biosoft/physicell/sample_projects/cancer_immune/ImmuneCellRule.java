package ru.biosoft.physicell.sample_projects.cancer_immune;

import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.PhysiCellUtilities;

public class ImmuneCellRule extends CustomCellRule
{
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        int attach_lifetime_i = pCell.custom_data.findVariableIndex( "attachment_lifetime" );

        if( phenotype.death.dead == true )
        {
            // the cell death functions don't automatically turn off custom functions, since those are part of mechanics. 
            // Let's just fully disable now. 
            pCell.functions.customCellRule = null;
            return;
        }

        // if I'm docked
        if( pCell.state.numberAttachedCells() > 0 )
        {
            // attempt to kill my attached cell
            Cell attached = pCell.state.attachedCells.iterator().next();///[0];
            boolean detach_me = false;

            if( attemptApoptosis( pCell, attached, dt ) )
            {
                triggerApoptosis( pCell, attached );
                detach_me = true;
            }

            // decide whether to detach 
            if( PhysiCellUtilities.UniformRandom() < dt / ( pCell.custom_data.get( attach_lifetime_i ) + 1e-15 ) )
                detach_me = true;

            // if I dettach, resume motile behavior 
            if( detach_me )
            {
                Cell.detachCells( pCell, attached );
                phenotype.motility.isMotile = true;
            }
            return;
        }

        // I'm not docked, look for cells nearby and try to docked if this returns non-NULL, we're now attached to a cell 
        if( checkNeighborsForAttachment( pCell, dt ) != null )
        {
            phenotype.motility.isMotile = false;
            return;
        }
        phenotype.motility.isMotile = true;
    }

    static boolean triggerApoptosis(Cell pAttacker, Cell pTarget)
    {
        if( pTarget.phenotype.death.dead )
            return false;
        pTarget.startDeath( pTarget.phenotype.death.findDeathModelIndex( "apoptosis" ) );
        return true;
    }

    static boolean attemptApoptosis(Cell pAttacker, Cell pTarget, double dt)
    {
        int oncoproteinIndex = pTarget.custom_data.findVariableIndex( "oncoprotein" );
        int killRateIndex = pAttacker.custom_data.findVariableIndex( "kill_rate" );

        double oncoproteinSaturation = pAttacker.custom_data.get( "oncoprotein_saturation" ); // 2.0; 
        double oncoproteinThreshold = pAttacker.custom_data.get( "oncoprotein_threshold" ); // 0.5; // 0.1; 
        double oncoproteinDifference = oncoproteinSaturation - oncoproteinThreshold;

        double targetOconoprotein = pTarget.custom_data.get( oncoproteinIndex );
        if( targetOconoprotein < oncoproteinThreshold )
            return false;

        double scale = ( targetOconoprotein - oncoproteinThreshold ) / oncoproteinDifference;
        scale = Math.min( scale, 1.0 );

        if( PhysiCellUtilities.UniformRandom() < pAttacker.custom_data.get( killRateIndex ) * scale * dt )
            return true;
        return false;
    }

    public static Cell checkNeighborsForAttachment(Cell pAttacker, double dt)
    {
        for( Cell nearbyCell : pAttacker.cells_in_my_container() )
        {
            if( nearbyCell != pAttacker )// don't try to kill yourself 
            {
                if( attemptAttachment( pAttacker, nearbyCell, dt ) )
                    return nearbyCell;
            }
        }
        return null;
    }

    static boolean attemptAttachment(Cell pAttacker, Cell pTarget, double dt)
    {
        double oncoprotein_saturation = pAttacker.custom_data.get( "oncoprotein_saturation" );
        double oncoprotein_threshold = pAttacker.custom_data.get( "oncoprotein_threshold" );
        double maxAttachmentDistance = pAttacker.custom_data.get( "max_attachment_distance" );
        double minAttachmentDistance = pAttacker.custom_data.get( "min_attachment_distance" );
        double targetOncoprotein = pTarget.custom_data.get( "oncoprotein" );
        if( targetOncoprotein > oncoprotein_threshold && !pTarget.phenotype.death.dead )
        {
            double distance = VectorUtil.dist( pTarget.position, pAttacker.position );
            if( distance > maxAttachmentDistance )
                return false;

            double attachRate = pAttacker.custom_data.get( "attachment_rate" );
            double scale = ( targetOncoprotein - oncoprotein_threshold ) / ( oncoprotein_saturation - oncoprotein_threshold );
            double distanceScale = ( maxAttachmentDistance - distance ) / ( maxAttachmentDistance - minAttachmentDistance );
            attachRate *= Math.min( scale, 1.0 ) * Math.min( distanceScale, 1.0 );
            if( PhysiCellUtilities.UniformRandom() < attachRate * dt )
                Cell.attachcCells( pAttacker, pTarget );
            return true;//TODO: should we return true only if attached successfully?
        }
        return false;
    }
}