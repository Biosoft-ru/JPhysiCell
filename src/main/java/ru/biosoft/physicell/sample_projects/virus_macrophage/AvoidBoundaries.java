package ru.biosoft.physicell.sample_projects.virus_macrophage;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;
import ru.biosoft.physicell.core.Phenotype;

/**
 *Adds velocity to steer clear of the boundaries 
 */
public class AvoidBoundaries extends CustomCellRule
{
    static double avoidZone = 25;
    static double avoidSpeed = -0.5; // must be negative 

    @Override
    public void execute(Cell cell, Phenotype phenotype, double dt)
    {
        Microenvironment m = cell.getMicroenvironment();
        double Xmin = m.mesh.boundingBox[0];
        double Ymin = m.mesh.boundingBox[1];
        double Zmin = m.mesh.boundingBox[2];
        double Xmax = m.mesh.boundingBox[3];
        double Ymax = m.mesh.boundingBox[4];
        double Zmax = m.mesh.boundingBox[5];

        boolean nearEdge = false;
        if( cell.position[0] < Xmin + avoidZone || cell.position[0] > Xmax - avoidZone )
            nearEdge = true;
        else if( cell.position[1] < Ymin + avoidZone || cell.position[1] > Ymax - avoidZone )
            nearEdge = true;
        else if( !m.options.simulate2D && ( cell.position[2] < Zmin + avoidZone || cell.position[2] > Zmax - avoidZone ) )
            //        {
            //            if( cell.position[2] < Zmin + avoidZone || cell.position[2] > Zmax - avoidZone )
                nearEdge = true;
        //        }

        if( nearEdge )
            //        {
            cell.velocity = VectorUtil.newProd( cell.position, avoidSpeed ); // move towards origin 
            //            VectorUtil.prod( cell.velocity, avoid_speed );
            //        }
    }
}