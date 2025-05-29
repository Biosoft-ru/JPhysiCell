package ru.biosoft.physicell.covid;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Model;

public class ExternalImmune
{
    double dC;
    double pT1;
    double pT2;
    double dT1;
    double dT2;
    double Tc0;
    double dDm;
    double sTh1;
    double pTh1;
    double dTh1;
    double mTh;
    double sTh2;
    double pTh2;
    double ro;
    double CD8_Tcell_recruitment_rate;
    double dB;
    double B0;
    double rB1;
    double h;
    double rB2;
    double pSc;
    double dS;
    double pAS;
    double dMc;
    ModelCovid model;

    private double x_min;
    private double x_max;
    private double y_min;
    private double y_max;

    private int nAb;
    private int nV;
    private Microenvironment microenvironment;

    public void setup(Model model)
    {
        this.model = (ModelCovid)model;
        dC = model.getParameterDouble( "TC_death_rate" );
        pT1 = model.getParameterDouble( "max_activation_TC" );
        pT2 = model.getParameterDouble( "half_max_activation_TC" );
        dT1 = model.getParameterDouble( "max_clearance_TC" );
        dT2 = model.getParameterDouble( "half_max_clearance_TC" );
        Tc0 = model.getParameterDouble( "TC_population_threshold" );
        dDm = model.getParameterDouble( "DM_decay" );
        sTh1 = model.getParameterDouble( "Th1_max_activation" );
        pTh1 = model.getParameterDouble( "Th1_damping" );
        dTh1 = model.getParameterDouble( "Th1_decay" );
        mTh = model.getParameterDouble( "Th_base_decay" );
        sTh2 = model.getParameterDouble( "Th2_self_feeback" );
        pTh2 = model.getParameterDouble( "Th2_max_conversion" );
        ro = model.getParameterDouble( "Th1_Th2_conversion_wieght" );
        CD8_Tcell_recruitment_rate = model.getParameterDouble( "T_Cell_Recruitment" );
        dB = model.getParameterDouble( "BCell_base_rate" );
        B0 = model.getParameterDouble( "BCell_base_value" );
        rB1 = model.getParameterDouble( "BCell_DC_proliferation" );
        h = model.getParameterDouble( "BCell_Th2_wieght_function" );
        rB2 = model.getParameterDouble( "BCell_damping" );
        pSc = model.getParameterDouble( "PCell_recuitment" );
        dS = model.getParameterDouble( "PCell_degradation" );
        pAS = model.getParameterDouble( "Ig_recuitment" );
        dMc = model.getParameterDouble( "Ig_degradation" );

        microenvironment = model.getMicroenvironment();
        x_min = microenvironment.mesh.boundingBox[0] + 1e-6;
        x_max = microenvironment.mesh.boundingBox[3] - 1e-6;
        y_min = microenvironment.mesh.boundingBox[1] + 1e-6;
        y_max = microenvironment.mesh.boundingBox[4] - 1e-6;

        nAb = microenvironment.findDensityIndex( "Ig" );
        nV = microenvironment.findDensityIndex( "virion" );
    }

    public void doStep(double dt)
    {


        double lypmh_scale = model.EPICOUNT / 500000;
        // actual model goes here 

        double[][] x = new double[4][9];//={{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}};//initialize x
        double[][] f = new double[4][9];//={{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}};//initialize f
        int j;

        // TC update
        double dR_TC = dC * Tc0;

        /*  // DM Tc recruitment
            double dR_TCD = pT1 * C[0]/immunevolume * C[1]/immunevolume / ( C[0]/immunevolume + pT2/immunevolume) ;
            
            // DM Tc decay
            double dR_TC16 = dT1 * C[0]/immunevolume * C[1]/immunevolume / ( C[0]/immunevolume + dT2/immunevolume) ;
            
            // TC decay
            double dR_TC14 = dC * C[1] / immunevolume ;
            
            // DM decay
            double dR_DM = dDm * C[0] / immunevolume; */



        //    extern std::vector<int>history;

        x[0][0] = ( model.DM + model.history.getBack() ) / lypmh_scale;
        x[0][1] = model.TC; //initial values
        x[0][2] = model.TH1; //initial values
        x[0][3] = model.TH2; //initial values
        x[0][4] = model.TCt / lypmh_scale;
        x[0][5] = model.Tht / lypmh_scale;
        x[0][6] = model.Bc;
        x[0][7] = model.Ps;
        x[0][8] = model.Ig / lypmh_scale;

        for( j = 0; j < 4; j++ )
        {
            f[j][0] = -dDm * x[j][0]; //DM
            f[j][1] = dR_TC - dC * x[j][1] + pT1 * ( ( 100000 - x[j][1] ) / ( 100000 ) ) * x[j][0] * x[j][1] / ( x[j][0] + pT2 )
                    - dT1 * x[j][0] * x[j][1] / ( x[j][0] + dT2 ); //Tc
            f[j][2] = ( sTh1 * x[j][2] ) / ( ( 1 + x[j][3] ) * ( 1 + x[j][3] ) )
                    + ( pTh1 * x[j][0] * x[j][2] * x[j][2] ) / ( ( 1 + x[j][3] ) * ( 1 + x[j][3] ) )
                    - ( dTh1 * x[j][0] * x[j][2] * x[j][2] * x[j][2] ) / ( 500 + x[j][3] ) - mTh * x[j][2]; //Th1
            f[j][3] = ( sTh2 * x[j][3] ) / ( 1 + x[j][3] )
                    + ( pTh2 * ( ro + x[j][2] ) * x[j][0] * x[j][3] * x[j][3] ) / ( ( 1 + x[j][3] ) * ( 1 + x[j][2] + x[j][3] ) )
                    - mTh * x[j][3]; //Th2
            f[j][4] = CD8_Tcell_recruitment_rate * x[j][1]; //CD8 export
            f[j][5] = 10 * CD8_Tcell_recruitment_rate * ( x[j][2] + x[j][3] ); //CD4 export
            f[j][6] = dB * B0 + rB1 * x[j][6] * ( x[j][0] + h * x[j][3] ) / ( x[j][0] + h * x[j][3] + rB2 ) - dB * x[j][6]
                    - 2 * pSc * x[j][6]; //B-Cell
            f[j][7] = pSc * x[j][6] - dS * x[j][7]; //P-Cell
            f[j][8] = pAS * x[j][7] - dMc * x[j][8]; //Ig
            if( j == 0 || j == 1 )
            {
                x[j + 1][0] = x[0][0] + dt / 2 * f[j][0]; //first and second x approximations
                x[j + 1][1] = x[0][1] + dt / 2 * f[j][1]; //first and second x approximations
                x[j + 1][2] = x[0][2] + dt / 2 * f[j][2]; //first and second x approximations
                x[j + 1][3] = x[0][3] + dt / 2 * f[j][3]; //first and second x approximations
                x[j + 1][4] = x[0][4] + dt / 2 * f[j][4]; //first and second x approximations
                x[j + 1][5] = x[0][5] + dt / 2 * f[j][5]; //first and second x approximations
                x[j + 1][6] = x[0][6] + dt / 2 * f[j][6]; //first and second x approximations
                x[j + 1][7] = x[0][7] + dt / 2 * f[j][7]; //first and second x approximations
                x[j + 1][8] = x[0][8] + dt / 2 * f[j][8]; //first and second x approximations
            }
            if( j == 2 )
            {
                x[j + 1][0] = x[0][0] + dt * f[j][0]; //third approximation
                x[j + 1][1] = x[0][1] + dt * f[j][1]; //third approximation
                x[j + 1][2] = x[0][2] + dt * f[j][2]; //third approximation
                x[j + 1][3] = x[0][3] + dt * f[j][3]; //third approximation
                x[j + 1][4] = x[0][4] + dt * f[j][4]; //third approximation
                x[j + 1][5] = x[0][5] + dt * f[j][5]; //third approximation
                x[j + 1][6] = x[0][6] + dt * f[j][6]; //third approximation
                x[j + 1][7] = x[0][7] + dt * f[j][7]; //third approximation
                x[j + 1][8] = x[0][8] + dt * f[j][8]; //third approximation
            }
        }

        model.DM = ( x[0][0] + dt * ( f[0][0] / 6 + f[1][0] / 3 + f[2][0] / 3 + f[3][0] / 6 ) ) * lypmh_scale;
        model.TC = x[0][1] + dt * ( f[0][1] / 6 + f[1][1] / 3 + f[2][1] / 3 + f[3][1] / 6 );
        model.TH1 = x[0][2] + dt * ( f[0][2] / 6 + f[1][2] / 3 + f[2][2] / 3 + f[3][2] / 6 ); //detirmine n+1
        model.TH2 = x[0][3] + dt * ( f[0][3] / 6 + f[1][3] / 3 + f[2][3] / 3 + f[3][3] / 6 ); //detirmine n+1
        model.TCt = ( x[0][4] + dt * ( f[0][4] / 6 + f[1][4] / 3 + f[2][4] / 3 + f[3][4] / 6 ) ) * lypmh_scale;
        model.Tht = ( x[0][5] + dt * ( f[0][5] / 6 + f[1][5] / 3 + f[2][5] / 3 + f[3][5] / 6 ) ) * lypmh_scale;
        model.Bc = x[0][6] + dt * ( f[0][6] / 6 + f[1][6] / 3 + f[2][6] / 3 + f[3][6] / 6 );
        model.Ps = x[0][7] + dt * ( f[0][7] / 6 + f[1][7] / 3 + f[2][7] / 3 + f[3][7] / 6 );
        model.Ig = ( x[0][8] + dt * ( f[0][8] / 6 + f[1][8] / 3 + f[2][8] / 3 + f[3][8] / 6 ) ) * lypmh_scale;

        //System.out.println(model.getCurrentTime()+" DM " + model.DM + " TC " + model.TC + " TH1 " + model.TH1 + " TH2 " + model.TH2 + " TCt " + model.TCt + " Tht " + model.Tht + " Bc " + model.Bc + " Ps " + model.Ps + " Ig " +model.Ig);
        
        double number_of_Ig = Math.floor( model.Ig );
        model.Ig -= number_of_Ig;

        if (model.Ig >=1 )
        System.out.println(  "Placing "+ number_of_Ig + " Ig ... " ); 
        if( number_of_Ig > 1000 )
        {
            number_of_Ig = 1000;
        }
        for( int n = 0; n < number_of_Ig; n++ )
        {
            // pick a random voxel 
            double[] position = new double[3];
            position[0] = x_min + ( x_max - x_min ) * model.getRNG().UniformRandom();
            position[1] = y_min + ( y_max - y_min ) * model.getRNG().UniformRandom();

            int m = model.getMicroenvironment().nearestVoxelIndex( position );

            microenvironment.get( m )[nAb] += 1.0 / microenvironment.mesh.dV;
        }
        //    #pragma omp parallel for
        for( int n = 0; n < microenvironment.numberVoxels(); n++ )
        {
            if( microenvironment.get( n )[nV] > 0 && microenvironment.get( n )[nAb] > 0 )
            {
                double rate = 1.5 * microenvironment.get( n )[nAb] * microenvironment.get( n )[nV] * dt; //rate is 1.5 after conversions for now - set to zero in no Ig cases
                if( rate < 0 )
                {
                    rate = 0;
                }
                if( rate > microenvironment.get( n )[nAb] )
                {
                    rate = microenvironment.get( n )[nAb];
                }
                if( rate > microenvironment.get( n )[nV] )
                {
                    rate = microenvironment.get( n )[nV];
                }
                microenvironment.get( n )[nAb] -= rate;
                microenvironment.get( n )[nV] -= rate;
            }
        }
    }
}
