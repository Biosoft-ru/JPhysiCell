package ru.biosoft.physicell.core;

import java.util.Set;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.standard.StandardModels;

public class CellContainerExperimental extends CellContainer
{
    public static String EXPERIMENTAL_CONTAINER_NAME = "Experimental";

    public CellContainerExperimental()
    {
        this.name = EXPERIMENTAL_CONTAINER_NAME;
    }

    @Override
    public void updateAllCells(Model model, double t, double phenotypeDT, double mechanicsDT, double diffusionDT) throws Exception
    {
        Microenvironment m = model.getMicroenvironment();
        double tTot = System.nanoTime();
        Set<Cell> cells = m.getAgents( Cell.class );

        //if it is the time for running cell cycle, do it!
        double timeSinceLastCycle = t - lastCellCycleTime;
        double phenotypeDTtolerance = 0.001 * phenotypeDT;
        double mechanicsDTtolerance = 0.001 * mechanicsDT;

        double time_since_last_mechanics = t - last_mechanics_time;

        cells.parallelStream().filter( c -> !c.isOutOfDomain ).forEach( cell -> {
            try
            {
                if( !cell.isOutOfDomain )
                {
                    cell.phenotype.secretion.advance( cell, cell.phenotype, diffusionDT );
                }

                if( !cell.isOutOfDomain && initialized )
                {
                    if( cell.phenotype.intracellular != null && cell.phenotype.intracellular.need_update( t ) )
                    {
                        if( ( cell.functions.pre_update_intracellular != null ) )
                            cell.functions.pre_update_intracellular.execute( cell, cell.phenotype, diffusionDT );

                        cell.phenotype.intracellular.update( cell, cell.phenotype, diffusionDT );

                        if( cell.functions.post_update_intracellular != null )
                            cell.functions.post_update_intracellular.execute( cell, cell.phenotype, diffusionDT );
                    }
                }
            }
            catch( Exception ex )
            {
                ex.printStackTrace();
            }

        } );


        if( !initialized )
        {
            time_since_last_mechanics = mechanicsDT;
            timeSinceLastCycle = phenotypeDT;
        }
        final double tslc = timeSinceLastCycle;

        if( Math.abs( timeSinceLastCycle - phenotypeDT ) < phenotypeDTtolerance )
        {
            cells.parallelStream().filter( c -> !c.isOutOfDomain ).forEach( cell -> {
                try
                {
                    cell.advanceBundledPhenotype( tslc, rulesEnabled );
                }
                catch( Exception ex )
                {

                }
            } );

            for( Cell cell : cellsReadyToDivide )
                cell.divide();

            for( Cell cell : cellsReadyToDie )
                cell.die();

            numDivisionsCurStep += cellsReadyToDivide.size();
            numDeathsCurStep += cellsReadyToDie.size();

            cellsReadyToDie.clear();
            cellsReadyToDivide.clear();
            lastCellCycleTime = t;
        }


        if( Math.abs( time_since_last_mechanics - mechanicsDT ) < mechanicsDTtolerance )
        {
            final double tslm = time_since_last_mechanics;

            if( m.options.calculate_gradients )
                m.computeAllGradientVectors();

            cells.parallelStream().forEach( cell -> {
                try
                {
                    if( !cell.isOutOfDomain )
                    {
                        if( cell.functions.contact != null )
                        {
                            for( Cell attached : cell.state.attachedCells )
                                cell.functions.contact.execute( cell, cell.phenotype, attached, attached.phenotype, tslm );
                        }

                        if( cell.functions.customCellRule != null )
                            cell.functions.customCellRule.execute( cell, cell.phenotype, tslm );

                        if( cell.functions.updateVelocity != null && cell.isMovable )
                            cell.functions.updateVelocity.execute( cell, cell.phenotype, tslm );
                    }
                    if( !model.disableAutomatedSpringAdhesions )
                    {
                        StandardModels.dynamic_spring_attachments( cell, cell.phenotype, tslm );
                        if( cell.isMovable )
                        {
                            for( Cell pC1 : cell.state.springAttachments )
                                StandardModels.standard_elastic_contact_function( cell, cell.phenotype, pC1, pC1.phenotype, tslm );
                        }
                    }
                    StandardModels.standard_cell_cell_interactions( cell, cell.phenotype, tslm );
                }
                catch( Exception ex )
                {
                    ex.printStackTrace();
                }
            } );


            for( Cell cell : cellsReadyToDie )
                cell.die();
            cellsReadyToDie.clear();

            cells.parallelStream().filter( c -> !c.isOutOfDomain && c.isMovable ).forEach( cell -> {
                cell.updatePosition( tslm );
                cell.updateVoxelInContainer();
            } );


            last_mechanics_time = t;
        }
        initialized = true;
        tTot = System.nanoTime() - tTot;
        tTotal += tTot;
    }
}
