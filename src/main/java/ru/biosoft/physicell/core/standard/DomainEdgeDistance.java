package ru.biosoft.physicell.core.standard;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.CellFunctions.DistanceCalculator;
import ru.biosoft.physicell.core.Phenotype;

public class DomainEdgeDistance extends DistanceCalculator
{
    double tolerance = 1e-7;
    double one_over_sqrt_2 = 0.70710678118;
    double one_over_sqrt_3 = 0.57735026919;

    public double execute(Cell pCell, Phenotype phenotype, double dummy)
    {
        Microenvironment m = pCell.getMicroenvironment();
        double min_distance = 9e99;
        int nearest_boundary = -1;

        // check against xL and xU
        double temp_distance = pCell.position[0] - m.mesh.boundingBox[0];
        if( temp_distance < min_distance )
        {
            min_distance = temp_distance;
            nearest_boundary = 0;
        }
        temp_distance = m.mesh.boundingBox[3] - pCell.position[0];
        if( temp_distance < min_distance )
        {
            min_distance = temp_distance;
            nearest_boundary = 1;
        }

        // check against yL and yU
        temp_distance = pCell.position[1] - m.mesh.boundingBox[1];
        if( temp_distance < min_distance )
        {
            min_distance = temp_distance;
            nearest_boundary = 2;
        }
        temp_distance = m.mesh.boundingBox[4] - pCell.position[1];
        if( temp_distance < min_distance )
        {
            min_distance = temp_distance;
            nearest_boundary = 3;
        }

        if( m.options.simulate2D == false )
        {
            // if in 3D, check against zL and zU
            temp_distance = pCell.position[2] - m.mesh.boundingBox[2];
            if( temp_distance < min_distance )
            {
                min_distance = temp_distance;
                nearest_boundary = 4;
            }
            temp_distance = m.mesh.boundingBox[5] - pCell.position[2];
            if( temp_distance < min_distance )
            {
                min_distance = temp_distance;
                nearest_boundary = 5;
            }

            // check for 3D exceptions 

            // lines 
            if( Math.abs( ( pCell.position[0] ) - ( pCell.position[1] ) ) < tolerance
                    && Math.abs( ( pCell.position[1] ) - ( pCell.position[2] ) ) < tolerance
                    && Math.abs( ( pCell.position[0] ) - ( pCell.position[2] ) ) < tolerance )
            {
                if( pCell.position[0] > 0 )
                {
                    if( pCell.position[0] > 0 && pCell.position[1] > 0 )
                    {
                        pCell.displacement = new double[] { -one_over_sqrt_3, -one_over_sqrt_3, -one_over_sqrt_3};
                    }
                    if( pCell.position[0] < 0 && pCell.position[1] > 0 )
                    {
                        pCell.displacement = new double[] {one_over_sqrt_3, -one_over_sqrt_3, -one_over_sqrt_3};
                    }

                    if( pCell.position[0] > 0 && pCell.position[1] < 0 )
                    {
                        pCell.displacement = new double[] { -one_over_sqrt_3, one_over_sqrt_3, -one_over_sqrt_3};
                    }
                    if( pCell.position[0] < 0 && pCell.position[1] < 0 )
                    {
                        pCell.displacement = new double[] {one_over_sqrt_3, one_over_sqrt_3, -one_over_sqrt_3};
                    }
                }
                else
                {
                    if( pCell.position[0] > 0 && pCell.position[1] > 0 )
                    {
                        pCell.displacement = new double[] { -one_over_sqrt_3, -one_over_sqrt_3, one_over_sqrt_3};
                    }
                    if( pCell.position[0] < 0 && pCell.position[1] > 0 )
                    {
                        pCell.displacement = new double[] {one_over_sqrt_3, -one_over_sqrt_3, one_over_sqrt_3};
                    }

                    if( pCell.position[0] > 0 && pCell.position[1] < 0 )
                    {
                        pCell.displacement = new double[] { -one_over_sqrt_3, one_over_sqrt_3, one_over_sqrt_3};
                    }
                    if( pCell.position[0] < 0 && pCell.position[1] < 0 )
                    {
                        pCell.displacement = new double[] {one_over_sqrt_3, one_over_sqrt_3, one_over_sqrt_3};
                    }
                }
                return min_distance;
            }

            // planes - let's not worry for today 

        }
        else
        {
            // check for 2D  exceptions 

            if( Math.abs( ( pCell.position[0] ) - ( pCell.position[1] ) ) < tolerance )
            {
                if( pCell.position[0] > 0 && pCell.position[1] > 0 )
                {
                    pCell.displacement = new double[] { -one_over_sqrt_2, -one_over_sqrt_2, 0};
                }
                if( pCell.position[0] < 0 && pCell.position[1] > 0 )
                {
                    pCell.displacement = new double[] {one_over_sqrt_2, -one_over_sqrt_2, 0};
                }

                if( pCell.position[0] > 0 && pCell.position[1] < 0 )
                {
                    pCell.displacement = new double[] { -one_over_sqrt_2, one_over_sqrt_2, 0};
                }
                if( pCell.position[0] < 0 && pCell.position[1] < 0 )
                {
                    pCell.displacement = new double[] {one_over_sqrt_2, one_over_sqrt_2, 0};
                }
                return min_distance;
            }
        }

        // no exceptions 
        switch( nearest_boundary )
        {
            case 0:
                pCell.displacement = new double[] {1, 0, 0};
                return min_distance;
            case 1:
                pCell.displacement = new double[] { -1, 0, 0};
                return min_distance;
            case 2:
                pCell.displacement = new double[] {0, 1, 0};
                return min_distance;
            case 3:
                pCell.displacement = new double[] {0, -1, 0};
                return min_distance;
            case 4:
                pCell.displacement = new double[] {0, 0, 1};
                return min_distance;
            case 5:
                pCell.displacement = new double[] {0, 0, -1};
                return min_distance;
            default:
                pCell.displacement = new double[] {0, 0, 0};
                return 9e99;
        }
        //            pCell.displacement = new double[] {0, 0, 0};
        //            return 9e99;
    }
}