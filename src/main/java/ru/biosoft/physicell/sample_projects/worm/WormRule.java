package ru.biosoft.physicell.sample_projects.worm;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.core.standard.Chemotaxis;

public class WormRule extends CustomCellRule
{
    private Model model;

    public WormRule(Model model)
    {
        this.model = model;
    }

    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        // bookkeeping 
        Microenvironment m = pCell.getMicroenvironment();
        int nSignal = m.findDensityIndex( "signal" );

        // look for cells to form attachments, if 0 attachments
        int number_of_attachments = pCell.state.numberAttachedCells();
        if( number_of_attachments == 0 )
        {
            for( Cell neighbor : pCell.nearby_interacting_cells() )
            {
                if( neighbor.state.numberAttachedCells() < neighbor.customData.get( "max_attachments" ) )
                {
                    Cell.attachcCells( neighbor, pCell );
                    number_of_attachments++;
                }
                if( number_of_attachments > pCell.customData.get( "max_attachments" ) )
                    break;
            }
        }

        if( number_of_attachments == 0 )
        {
            pCell.functions.updateMigration = new Chemotaxis();
        }
        else if( number_of_attachments == 1 ) // if 1 attachment, do some logic  
        {
            // constant expression in end cells 
            pCell.customData.set( "head", pCell.customData.get( "head_initial" ) );

            // am I the head? 
            boolean head = false;
            if( pCell.customData.get( "head" ) > pCell.state.attachedCells.iterator().next().customData.get( "head" ) )
                head = true;

            //            if( head )
            pCell.functions.updateMigration = head ? new HeadMigration( model ) : new TailMigration( model );
            //            else
            //                pCell.functions.updateMigration = new TailMigration( model );

            phenotype.secretion.secretionRates[nSignal] = 100;
        }
        else if( number_of_attachments > 1 ) // if 2 or more attachments, use middle 
        {
            pCell.functions.updateMigration = new MiddleMigration( model );
            phenotype.secretion.secretionRates[nSignal] = 1;
        }
    }
}