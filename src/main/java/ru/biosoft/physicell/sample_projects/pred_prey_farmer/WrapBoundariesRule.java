package ru.biosoft.physicell.sample_projects.pred_prey_farmer;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;
import ru.biosoft.physicell.core.Phenotype;

public class WrapBoundariesRule extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        wrapBoundaries( pCell );
    }

    public static void wrapBoundaries(Cell pCell)
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        // add velocity to steer clear of the boundaries 
        double xMin = microenvironment.mesh.boundingBox[0];
        double yMin = microenvironment.mesh.boundingBox[1];
        double zMin = microenvironment.mesh.boundingBox[2];

        double xMax = microenvironment.mesh.boundingBox[3];
        double yMax = microenvironment.mesh.boundingBox[4];
        double zMax = microenvironment.mesh.boundingBox[5];

        double avoidZone = 20;

        boolean setupDone = false;
        if( setupDone == false )
        {
            xMax -= avoidZone;
            xMin += avoidZone;

            yMax -= avoidZone;
            yMin += avoidZone;

            setupDone = true;
        }

        boolean wrapped = false;

        double[] p = pCell.position.clone();
        double delta;


        while( p[0] < xMin )
        {
            delta = xMin - p[0];
            p[0] = xMax - delta;
            wrapped = true;
        }
        while( p[0] > xMax )
        {
            delta = p[0] - xMax;
            p[0] = xMin + delta;
            wrapped = true;
        }

        while( p[1] < yMin )
        {
            delta = yMin - p[1];
            p[1] = yMax - delta;
            wrapped = true;
        }
        while( p[1] > yMax )
        {
            delta = p[1] - yMax;
            p[1] = yMin + delta;
            wrapped = true;
        }

        if( !microenvironment.options.simulate2D )
        {
            while( p[2] < zMin )
            {
                delta = zMin - p[2];
                p[2] = zMax - delta;
                wrapped = true;
            }
            while( p[2] > zMax )
            {
                delta = p[2] - zMax;
                p[2] = zMin + delta;
                wrapped = true;
            }
        }

        if( !wrapped )
        {
            pCell.assignPosition( p );
        }
    }

    @Override
    public String display()
    {
        return "Wrap microenvironment boundaries";
    }
}