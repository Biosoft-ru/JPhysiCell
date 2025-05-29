package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Model;

public class ReceptorDynamics
{
    int nV_external;
    int nA_external;

    int nV_internal;
    int nA_internal;

    // bookkeeping -- find custom data we need 

    int nR_EU;
    int nR_EB;
    int nR_IU;
    int nR_IB;

    int nR_bind;
    int nR_endo;
    int nR_release;
    int nR_recycle;

    int lung_epithelial_type;
    int ignore_smoothing_flag;

    double x_min;
    double x_max;
    double y_min;
    double y_max;

    private Microenvironment m;
    private Model model;

    public void setup(Model model)
    {
        this.model = model;
        m = model.getMicroenvironment();
        ignore_smoothing_flag = model.getParameterInt( "ignore_smoothing_flag" );

        nV_external = m.findDensityIndex( "virion" );
        nA_external = m.findDensityIndex( "assembled virion" );

        lung_epithelial_type = model.getCellDefinition( "lung epithelium" ).type;

        x_min = m.mesh.boundingBox[0];
        x_max = m.mesh.boundingBox[3];
        y_min = m.mesh.boundingBox[1];
        y_max = m.mesh.boundingBox[4];
        // submodel_registry.register_model( receptor_dynamics_info ); 
        //        receptor_dynamics_info.register_model();   
    }

    void receptor_dynamics_model(Cell pCell, double dt)
    {
        double before = pCell.nearest_density_vector()[nV_external];
        double add1 = 0;
        double add2 = 0;
        double sub1 = 0;
        double sub2 = 0;
        double sub3 = 0;
        nV_internal = pCell.customData.findVariableIndex( "virion" );
        nA_internal = pCell.customData.findVariableIndex( "assembled_virion" );

        nR_EU = pCell.customData.findVariableIndex( "unbound_external_ACE2" );
        nR_EB = pCell.customData.findVariableIndex( "bound_external_ACE2" );
        nR_IU = pCell.customData.findVariableIndex( "unbound_internal_ACE2" );
        nR_IB = pCell.customData.findVariableIndex( "bound_internal_ACE2" );

        nR_bind = pCell.customData.findVariableIndex( "ACE2_binding_rate" );
        nR_endo = pCell.customData.findVariableIndex( "ACE2_endocytosis_rate" );
        nR_release = pCell.customData.findVariableIndex( "ACE2_cargo_release_rate" );
        nR_recycle = pCell.customData.findVariableIndex( "ACE2_recycling_rate" );

        // do nothing if dead 
        if( pCell.phenotype.death.dead )
            return;

        // if not lung epithelium, do nothing 
        if( pCell.type != lung_epithelial_type )
            return;

        // actual model goes here 
        double[][] x = new double[4][6];
        double[][] f = new double[4][6];
        int j;//initialize counter

        //initial values for RK4
        x[0][0] = pCell.customData.get( nR_EU );
        x[0][1] = pCell.customData.get( nR_EB );
        x[0][2] = pCell.customData.get( nR_IB );
        x[0][3] = pCell.customData.get( nR_IU );
        x[0][4] = pCell.customData.get( nV_internal );
        x[0][5] = 0;

        /*  // internalize
        double dR_IB = pCell.custom_data[nR_endo]*pCell.custom_data[nR_EB];   
        // viral release from endosomes     
        double dR_IU = pCell.custom_data[nR_release]*pCell.custom_data[nR_IB];    
        // receptor recycling   
        double dR_EU = pCell.custom_data[nR_recycle]*pCell.custom_data[nR_IU];
         */
        //int ignore_smoothing_flag=1;

        double dt_bind = dt * pCell.customData.get( nR_bind ) * pCell.nearest_density_vector()[nV_external] * pCell.phenotype.volume.total
                * pCell.customData.get( nR_EU ); //use FE to find what loop to enter

        if( dt_bind < 1 )
        {
            //SOLVE ODE BEFORE STOCHASTIC PORTION
            for( j = 0; j < 4; j++ )
            {
                f[j][0] = pCell.customData.get( nR_recycle ) * x[j][3]; //define SPECIAL function
                f[j][1] = -pCell.customData.get( nR_endo ) * x[j][1] - 0 * x[j][1]; //define SPECIAL function
                f[j][2] = pCell.customData.get( nR_endo ) * x[j][1] - pCell.customData.get( nR_release ) * x[j][2]; //define function
                f[j][3] = pCell.customData.get( nR_release ) * x[j][2] - pCell.customData.get( nR_recycle ) * x[j][3]; //define function
                f[j][4] = pCell.customData.get( nR_release ) * x[j][2]; //define function
                f[j][5] = 0 * x[j][1]; //counter for export
                if( j == 0 || j == 1 )
                {
                    x[j + 1][0] = x[0][0] + dt / 2 * f[j][0]; //first and second x approximations
                    x[j + 1][1] = x[0][1] + dt / 2 * f[j][1]; //first and second x approximations
                    x[j + 1][2] = x[0][2] + dt / 2 * f[j][2]; //first and second x approximations
                    x[j + 1][3] = x[0][3] + dt / 2 * f[j][3]; //first and second x approximations
                    x[j + 1][4] = x[0][4] + dt / 2 * f[j][4]; //first and second x approximations
                    x[j + 1][5] = x[0][5] + dt / 2 * f[j][5]; //first and second x approximations
                }
                if( j == 2 )
                {
                    x[j + 1][0] = x[0][0] + dt * f[j][0]; //third approximation
                    x[j + 1][1] = x[0][1] + dt * f[j][1]; //third approximation
                    x[j + 1][2] = x[0][2] + dt * f[j][2]; //third approximation
                    x[j + 1][3] = x[0][3] + dt * f[j][3]; //third approximation
                    x[j + 1][4] = x[0][4] + dt * f[j][4]; //third approximation
                    x[j + 1][5] = x[0][5] + dt * f[j][5]; //third approximation
                }
            }

            pCell.customData.set( nR_EU, x[0][0] + dt * ( f[0][0] / 6 + f[1][0] / 3 + f[2][0] / 3 + f[3][0] / 6 ) );
            pCell.customData.set( nR_EB, x[0][1] + dt * ( f[0][1] / 6 + f[1][1] / 3 + f[2][1] / 3 + f[3][1] / 6 ) );
            pCell.customData.set( nR_IB, x[0][2] + dt * ( f[0][2] / 6 + f[1][2] / 3 + f[2][2] / 3 + f[3][2] / 6 ) ); //detirmine n+1
            pCell.customData.set( nR_IU, x[0][3] + dt * ( f[0][3] / 6 + f[1][3] / 3 + f[2][3] / 3 + f[3][3] / 6 ) ); //detirmine n+1
            pCell.customData.set( nV_internal, x[0][4] + dt * ( f[0][4] / 6 + f[1][4] / 3 + f[2][4] / 3 + f[3][4] / 6 ) ); //detirmine n+1

            add1 = dt * ( f[0][5] / 6 + f[1][5] / 3 + f[2][5] / 3 + f[3][5] / 6 ) / m.mesh.dV;
            pCell.nearest_density_vector()[nV_external] += dt * ( f[0][5] / 6 + f[1][5] / 3 + f[2][5] / 3 + f[3][5] / 6 ) / m.mesh.dV;

            //START STOCHASTIC PORTION
            if( dt_bind > 0 && model.getRNG().UniformRandom() <= dt_bind )
            {
                // don't attach more virus than the available unbound receptor
                if( pCell.customData.get( nR_EU ) >= 1 )
                {
                    pCell.customData.modify( nR_EU, -1 );
                    pCell.customData.modify( nR_EB, 1 );

                   
                    sub1 = -1.0/ m.mesh.dV;
                    if( ignore_smoothing_flag > 0.5 )
                    {
                        pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                    }
                    else
                    {
                        if( pCell.nearest_density_vector()[nV_external] >= 0.5 / m.mesh.dV )
                        {
                            pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                        }
                        else
                        {
                            double p0 = pCell.position[0];
                            double p1 = pCell.position[1];
                            double p2 = pCell.position[2];

                            double[] dummypos = new double[] {p0 + 20, p1, p2};
                            double[] dummypos0 = new double[] {p0 - 20, p1, p2};
                            double[] dummypos1 = new double[] {p0, p1 + 20, p2};
                            double[] dummypos2 = new double[] {p0, p1 - 20, p2};
                            double[] dummypos00 = new double[] {p0 - 20, p1 + 20, p2};
                            double[] dummypos01 = new double[] {p0 + 20, p1 + 20, p2};
                            double[] dummypos11 = new double[] {p0 + 20, p1 - 20, p2};
                            double[] dummypos10 = new double[] {p0 - 20, p1 - 20, p2};

                            if( dummypos[0] > x_max )
                            {
                                if( dummypos1[1] > y_max )
                                {
                                    pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                                }
                                else if( dummypos2[1] < y_min )
                                {
                                    pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                                }
                                else
                                {
                                    pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                                }
                            }
                            else if( dummypos0[0] < x_min )
                            {
                                if( dummypos1[1] > y_max )
                                {
                                    pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                                }
                                else if( dummypos2[1] < y_min )
                                {
                                    pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                                }
                                else
                                {
                                    pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                                }
                            }
                            else
                            {
                                if( dummypos1[1] > y_max )
                                {
                                    pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                                }
                                else if( dummypos2[1] < y_min )
                                {
                                    pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                                }
                                else
                                {
                                    pCell.nearest_density_vector()[nV_external] -= 0.25 / m.mesh.dV;
                                    sub1 = -0.25/ m.mesh.dV;
                                    m.nearestDensity( m.nearestVoxelIndex( dummypos ) )[nV_external] -= 0.125 / m.mesh.dV;
                                    m.nearestDensity( m.nearestVoxelIndex( dummypos0 ) )[nV_external] -= 0.125 / m.mesh.dV;
                                    m.nearestDensity( m.nearestVoxelIndex( dummypos1 ) )[nV_external] -= 0.125 / m.mesh.dV;
                                    m.nearestDensity( m.nearestVoxelIndex( dummypos2 ) )[nV_external] -= 0.125 / m.mesh.dV;
                                    m.nearestDensity( m.nearestVoxelIndex( dummypos00 ) )[nV_external] -= 0.0625 / m.mesh.dV;
                                    m.nearestDensity( m.nearestVoxelIndex( dummypos01 ) )[nV_external] -= 0.0625 / m.mesh.dV;
                                    m.nearestDensity( m.nearestVoxelIndex( dummypos11 ) )[nV_external] -= 0.0625 / m.mesh.dV;
                                    m.nearestDensity( m.nearestVoxelIndex( dummypos10 ) )[nV_external] -= 0.0625 / m.mesh.dV;
                                }
                            }
                        }
                    }
                }
            }

        }

        if( dt_bind > 1 )
        {
            //SOLVE ODE BEFORE STOCHASTIC PORTION
            //THIS RK4 METHOD SOLVES THE BINDING IN A DETERMINISTIC FASHON

            for( j = 0; j < 4; j++ )
            {
                f[j][0] = pCell.customData.get( nR_recycle ) * x[j][3]; //define SPECIAL function
                f[j][1] = -pCell.customData.get( nR_endo ) * x[j][1] - 0 * x[j][1]; //define SPECIAL function
                f[j][2] = pCell.customData.get( nR_endo ) * x[j][1] - pCell.customData.get( nR_release ) * x[j][2]; //define function
                f[j][3] = pCell.customData.get( nR_release ) * x[j][2] - pCell.customData.get( nR_recycle ) * x[j][3]; //define function
                f[j][4] = pCell.customData.get( nR_release ) * x[j][2]; //define function
                f[j][5] = 0 * x[j][1]; //counter for export
                if( j == 0 || j == 1 )
                {
                    x[j + 1][0] = x[0][0] + dt / 2 * f[j][0]; //first and second x approximations
                    x[j + 1][1] = x[0][1] + dt / 2 * f[j][1]; //first and second x approximations
                    x[j + 1][2] = x[0][2] + dt / 2 * f[j][2]; //first and second x approximations
                    x[j + 1][3] = x[0][3] + dt / 2 * f[j][3]; //first and second x approximations
                    x[j + 1][4] = x[0][4] + dt / 2 * f[j][4]; //first and second x approximations
                    x[j + 1][5] = x[0][5] + dt / 2 * f[j][5]; //first and second x approximations
                }
                if( j == 2 )
                {
                    x[j + 1][0] = x[0][0] + dt * f[j][0]; //third approximation
                    x[j + 1][1] = x[0][1] + dt * f[j][1]; //third approximation
                    x[j + 1][2] = x[0][2] + dt * f[j][2]; //third approximation
                    x[j + 1][3] = x[0][3] + dt * f[j][3]; //third approximation
                    x[j + 1][4] = x[0][4] + dt * f[j][4]; //third approximation
                    x[j + 1][5] = x[0][5] + dt * f[j][5]; //third approximation
                }
            }

            pCell.customData.set( nR_EU, x[0][0] + dt * ( f[0][0] / 6 + f[1][0] / 3 + f[2][0] / 3 + f[3][0] / 6 ) );
            pCell.customData.set( nR_EB, x[0][1] + dt * ( f[0][1] / 6 + f[1][1] / 3 + f[2][1] / 3 + f[3][1] / 6 ) );
            pCell.customData.set( nR_IB, x[0][2] + dt * ( f[0][2] / 6 + f[1][2] / 3 + f[2][2] / 3 + f[3][2] / 6 ) ); //detirmine n+1
            pCell.customData.set( nR_IU, x[0][3] + dt * ( f[0][3] / 6 + f[1][3] / 3 + f[2][3] / 3 + f[3][3] / 6 ) ); //detirmine n+1
            pCell.customData.set( nV_internal, x[0][4] + dt * ( f[0][4] / 6 + f[1][4] / 3 + f[2][4] / 3 + f[3][4] / 6 ) ); //detirmine n+1

            //attempting proper integration to stochastic portion, still needs some thought
            double alpha =  dt_bind;
            double n_virion = pCell.nearest_density_vector()[nV_external] * m.mesh.dV;

            pCell.nearest_density_vector()[nV_external] += dt * ( f[0][5] / 6 + f[1][5] / 3 + f[2][5] / 3 + f[3][5] / 6 ) / m.mesh.dV;
            add2 = dt * ( f[0][5] / 6 + f[1][5] / 3 + f[2][5] / 3 + f[3][5] / 6 ) / m.mesh.dV;
            //limit to number of virons in a voxel
            if( alpha > n_virion )
            {
                alpha = n_virion;
            }
            double alpha1 = Math.floor( alpha );
            double alpha2 = alpha - alpha1;

            //limit to number of unbound receptor on cell surface
            if( alpha1 > pCell.customData.get( nR_EU ) )
            {
                alpha1 = pCell.customData.get( nR_EU );
            }

            //STOCHASTIC PORTION
            pCell.customData.modify( nR_EU, -alpha1 );
            pCell.customData.modify( nR_EB, alpha1 );

            pCell.nearest_density_vector()[nV_external] -= alpha1 / m.mesh.dV;

            sub2 = -alpha1 / m.mesh.dV;
            if( model.getRNG().UniformRandom() <= alpha2 )
            {
                // don't attach more virus than the available unbound receptor
                if( pCell.customData.get( nR_EU ) >= 1 )
                {
                    pCell.customData.modify( nR_EU, -1 );
                    pCell.customData.modify( nR_EB, 1 );
                    pCell.nearest_density_vector()[nV_external] -= 1.0 / m.mesh.dV;
                    sub3 = -1.0 / m.mesh.dV;
                }

            }
        }
        if( pCell.nearest_density_vector()[nV_external] < 0 )
        {
//            System.out.println( "" );
        }
    }

    public void doStep(double dt)
    {
        for( Cell cell : model.getMicroenvironment().getAgents( Cell.class ) )
        {
            if( !cell.phenotype.death.dead )
                receptor_dynamics_model( cell, dt );
        }
    }
}
