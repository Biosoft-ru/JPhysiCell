package ru.biosoft.physicell.micromets;

import ru.biosoft.physicell.micromets.MicrometsModel.LymphNode;

public class external_immue
{
    MicrometsModel model;
    LymphNode lympNode;
    History history;

    public external_immue(MicrometsModel model)
    {
        this.model = model;
        this.lympNode = model.lympNode;
        this.history = model.history;
    }
    
    void external_immune_model(double dt)
    {
        double dC = model.getParameterDouble( "TC_death_rate" );
        double pT1 = model.getParameterDouble( "activation_rate_TC" );
        double pT2 = model.getParameterDouble( "half_max_activation_TC" );
        double dT1 = model.getParameterDouble( "clearance_rate_TC" );
        double dT2 = model.getParameterDouble( "half_max_clearance_TC" );
        double Tc0 = model.getParameterDouble( "TC_population_threshold" );
        double dDm = model.getParameterDouble( "DM_decay" );
        double sTh1 = model.getParameterDouble( "Th1_max_activation" );
        double pTh1 = model.getParameterDouble( "Th1_damping" ); //8.3e-6;
        double dTh1 = model.getParameterDouble( "Th1_decay" );
        double mTh = model.getParameterDouble( "Th_base_decay" ); //1.5e-5;
        double sTh2 = model.getParameterDouble( "Th2_self_feeback" ); //3e-5
        double pTh2 = model.getParameterDouble( "Th2_max_conversion" );
        double ro = model.getParameterDouble( "Th1_Th2_conversion_weight" );
        double CD8_Tcell_recruitment_rate = model.getParameterDouble( "T_Cell_Recruitment" );

        double lypmh_scale = model.lympNode.GridCOUNT / 5e5;

        // actual model goes here
        double[][] x = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};//initialize x
        double[][] f = {{0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}};//initialize f
        int j;

        // TC update
        double dR_TC = dC * Tc0;

        //        extern std::vector<int>history;

        x[0][0] = ( lympNode.DM + history.getBack() ) / lypmh_scale;
        x[0][1] = lympNode.TC; //initial values
        x[0][2] = lympNode.TH1; //initial values
        x[0][3] = lympNode.TH2; //initial values
        x[0][4] = lympNode.TCt / lypmh_scale;
        x[0][5] = lympNode.Tht / lypmh_scale;

        for( j = 0; j < 4; j++ )
        {
            f[j][0] = -dDm * x[j][0]; //define function
            f[j][1] = dR_TC - dC * x[j][1] + pT1 * ( ( 1000000 - x[j][1] ) / ( 1000000 ) ) * x[j][0] * x[j][1] / ( x[j][0] + pT2 )
                    - dT1 * x[j][0] * x[j][1] / ( x[j][0] + dT2 );
            f[j][2] = ( sTh1 * x[j][2] ) / ( ( 1 + x[j][3] ) * ( 1 + x[j][3] ) )
                    + ( pTh1 * x[j][0] * x[j][2] * x[j][2] ) / ( ( 1 + x[j][3] ) * ( 1 + x[j][3] ) )
                    - ( dTh1 * x[j][0] * x[j][2] * x[j][2] * x[j][2] ) / ( 500 + x[j][3] ) - mTh * x[j][2]; //define function
            f[j][3] = ( sTh2 * x[j][3] ) / ( 1 + x[j][3] )
                    + ( pTh2 * ( ro + x[j][2] ) * x[j][0] * x[j][3] * x[j][3] ) / ( ( 1 + x[j][3] ) * ( 1 + x[j][2] + x[j][3] ) )
                    - mTh * x[j][3]; //define function
            f[j][4] = CD8_Tcell_recruitment_rate * x[j][1]; //define function
            f[j][5] = CD8_Tcell_recruitment_rate * ( x[j][2] + x[j][3] ); //define function
            
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

        lympNode.DM = ( x[0][0] + dt * ( f[0][0] / 6 + f[1][0] / 3 + f[2][0] / 3 + f[3][0] / 6 ) ) * lypmh_scale;
        lympNode.TC = x[0][1] + dt * ( f[0][1] / 6 + f[1][1] / 3 + f[2][1] / 3 + f[3][1] / 6 );
        lympNode.TH1 = x[0][2] + dt * ( f[0][2] / 6 + f[1][2] / 3 + f[2][2] / 3 + f[3][2] / 6 ); //detirmine n+1
        lympNode.TH2 = x[0][3] + dt * ( f[0][3] / 6 + f[1][3] / 3 + f[2][3] / 3 + f[3][3] / 6 ); //detirmine n+1
        lympNode.TCt = ( x[0][4] + dt * ( f[0][4] / 6 + f[1][4] / 3 + f[2][4] / 3 + f[3][4] / 6 ) ) * lypmh_scale;
        lympNode.Tht = ( x[0][5] + dt * ( f[0][5] / 6 + f[1][5] / 3 + f[2][5] / 3 + f[3][5] / 6 ) ) * lypmh_scale;

        return;
    }
}
