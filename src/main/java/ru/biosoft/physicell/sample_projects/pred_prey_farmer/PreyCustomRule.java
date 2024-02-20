package ru.biosoft.physicell.sample_projects.pred_prey_farmer;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.biofvm.VectorUtil;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.custom_cell_rule;
import ru.biosoft.physicell.core.Phenotype;

public class PreyCustomRule implements custom_cell_rule
{
    @Override
    public void execute(Cell pCell, Phenotype phenotype, double dt)
    {
        wrap_boundaries( pCell );
    }

    public static void wrap_boundaries(Cell pCell)
    {

        avoid_boundaries( pCell );
        return;

        //    Microenvironment microenvironment = pCell.getMicroenvironment();
        //  // add velocity to steer clear of the boundaries 
        //  double Xmin = microenvironment.mesh.boundingBox[0]; 
        //  double Ymin = microenvironment.mesh.boundingBox[1]; 
        //  double Zmin = microenvironment.mesh.boundingBox[2]; 
        //
        //  double Xmax = microenvironment.mesh.boundingBox[3]; 
        //  double Ymax = microenvironment.mesh.boundingBox[4]; 
        //  double Zmax = microenvironment.mesh.boundingBox[5]; 
        //  
        //  double avoid_zone = 20; 
        //
        //  boolean setup_done = false; 
        //  if( setup_done == false )
        //  {
        //      Xmax -= avoid_zone; 
        //      Xmin += avoid_zone; 
        //      
        //      Ymax -= avoid_zone;
        //      Ymin += avoid_zone; 
        //      
        //      setup_done = true; 
        //  }
        //  
        //  boolean wrapped = false; 
        //  
        //  double[] p = pCell.position;
        //  double Delta;
        //
        //
        //  while( p[0] < Xmin )
        //  {
        //      Delta = Xmin - p[0]; 
        //      p[0] = Xmax - Delta; 
        //      wrapped = true; 
        //  }
        //  while( p[0] > Xmax )
        //  {
        //      Delta = p[0] - Xmax; 
        //      p[0] = Xmin + Delta; 
        //      wrapped = true; 
        //  }
        //  
        //  while( p[1] < Ymin )
        //  {
        //      Delta = Ymin - p[1]; 
        //      p[1] = Ymax - Delta; 
        //      wrapped = true; 
        //  }
        //  while( p[1] > Ymax )
        //  {
        //      Delta = p[1] - Ymax; 
        //      p[1] = Ymin + Delta; 
        //      wrapped = true; 
        //  }
        //
        //    if( microenvironment.options.simulate2D == false )
        //  {
        //      while( p[2] < Zmin )
        //      {
        //          Delta = Zmin - p[2]; 
        //          p[2] = Zmax - Delta; 
        //          wrapped = true; 
        //      }
        //      while( p[2] > Zmax )
        //      {
        //          Delta = p[2] - Zmax; 
        //          p[2] = Zmin + Delta; 
        //          wrapped = true; 
        //      }
        //  }
        //
        //  if( wrapped == true ) 
        //  { 
        //        //        #pragma omp critical 
        //        {
        //            pCell.assignPosition( p );
        //        }
        //  } 
    }

    public static void avoid_boundaries(Cell pCell)
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
            //        pCell.velocity *= avoid_speed; // move towards origin 
        }
    }
}