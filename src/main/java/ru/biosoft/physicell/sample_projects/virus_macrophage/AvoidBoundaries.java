package ru.biosoft.physicell.sample_projects.virus_macrophage;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.CustomCellRule;
import ru.biosoft.physicell.core.Phenotype;

public class AvoidBoundaries extends CustomCellRule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        Microenvironment microenvironment = pCell.getMicroenvironment();
        // add velocity to steer clear of the boundaries 
        double Xmin = microenvironment.mesh.boundingBox[0];
        double Ymin = microenvironment.mesh.boundingBox[1];
        double Zmin = microenvironment.mesh.boundingBox[2];

        double Xmax = microenvironment.mesh.boundingBox[3];
        double Ymax = microenvironment.mesh.boundingBox[4];
        double Zmax = microenvironment.mesh.boundingBox[5];

        double avoid_zone = 25;
        double avoid_speed = -0.5; // must be negative 

        // near edge: 
        boolean near_edge = false;
        if( pCell.position[0] < Xmin + avoid_zone || pCell.position[0] > Xmax - avoid_zone )
        {
            near_edge = true;
        }

        if( pCell.position[1] < Ymin + avoid_zone || pCell.position[1] > Ymax - avoid_zone )
        {
            near_edge = true;
        }

        if( microenvironment.options.simulate2D == false )
        {
            if( pCell.position[2] < Zmin + avoid_zone || pCell.position[2] > Zmax - avoid_zone )
            {
                near_edge = true;
            }
        }

        if( near_edge )
        {
            pCell.velocity = pCell.position; // move towards origin 
            VectorUtil.prod( pCell.velocity, avoid_speed );
            //            pCell.velocity *= avoid_speed; // move towards origin 
        }
    }
}