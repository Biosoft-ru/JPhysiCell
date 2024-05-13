package ru.biosoft.physicell.ode;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Intracellular;
import ru.biosoft.physicell.core.Phenotype;

public abstract class IntracellularODE extends Intracellular
{
    protected Map<String, String> phenotype_species = new HashMap<>();

    @Override
    public int updatePhenotypeParameters(Microenvironment m, Phenotype phenotype) throws Exception
    {
        for( Entry<String, String> elm : phenotype_species.entrySet() )
        {
            String first = elm.getKey();
            // motility params
            if( first.startsWith( "m" ) )
            {
                if( first.equals( "mms" ) )
                {
                    phenotype.motility.migrationSpeed = phenotype.intracellular.getParameterValue( elm.getValue() );
                }
                else if( first.equals( "mpt" ) )
                {
                    phenotype.motility.persistenceTime = phenotype.intracellular.getParameterValue( elm.getValue() );
                }
                else if( first.equals( "mmb" ) )
                {
                    phenotype.motility.migrationBias = phenotype.intracellular.getParameterValue( elm.getValue() );
                }
                else
                {
                }
            }
            // death params
            else if( first.startsWith( "d" ) )
            {
                if( first.equals( "da" ) )
                {
                    phenotype.death.rates.set( 0, phenotype.intracellular.getParameterValue( elm.getValue() ) );
                }
                else if( first.equals( "dn" ) )
                {
                    phenotype.death.rates.set( 1, phenotype.intracellular.getParameterValue( elm.getValue() ) );
                }
                else
                {
                }
            }
            // secretion params
            else if( first.startsWith( "s" ) )
            {
                // parsing attribute and getting substrate name
                //                String s = first;
                //                String delimiter = "_";

                //                int pos = 0;
                //                String token;
                //                while( ( pos = s.find( delimiter ) ) != String::npos )
                //                {
                //                    token = s.substr( 0, pos );
                //                    s.erase( 0, pos + delimiter.length() );
                //                }
                String[] tokens = first.split( "_" );
                int sub_index = m.findDensityIndex( tokens[1] );

                //transport types
                //uptake rate
                if( first.substring( 0, 3 ).equals( "sur" ) )
                {
                    //std::cout + sub_index + \n;
                    //std::cout + "Before sur1 : " + phenotype.secretion.uptake_rates[sub_index] + \n;
                    phenotype.secretion.uptakeRates[1] = phenotype.intracellular.getParameterValue( elm.getValue() );
                    //std::cout + "After sur1 : " + phenotype.secretion.uptake_rates[sub_index] + \n;
                }
                //secretion rate
                else if( first.substring( 0, 3 ).equals( "ssr" ) )
                {
                    phenotype.secretion.secretionRates[sub_index] = phenotype.intracellular.getParameterValue( elm.getValue() );
                }
                //secretion density
                else if( first.substring( 0, 3 ).equals( "ssd" ) )
                {
                    phenotype.secretion.saturationDensities[sub_index] = phenotype.intracellular.getParameterValue( elm.getValue() );
                }
                //net export rate
                else if( first.substring( 0, 3 ).equals( "ser" ) )
                {
                    phenotype.secretion.netExportRates[sub_index] = phenotype.intracellular.getParameterValue( elm.getValue() );
                }
                else
                {
                }
            }

            // cycle params
            else if( first.startsWith( "c" ) )
            {
                if( first.substring( 0, 3 ).equals( "ctr" ) )
                {
                    // parsing attribute and getting substrate name
                    //                    String s = first;
                    String delimiter = "_";

                    //                    int pos = 0;
                    //                    String token;
                    //                    int counter = 0;
                    //                    int start_index;
                    String[] tokens = first.split( delimiter );
                    //                    while( ( pos = s.find( delimiter ) ) != String::npos )
                    //                    {
                    //                        token = s.substring( 0, pos );
                    //                        //std::cout + counter + " : "+ token + \n;
                    //                        if( counter == 1 )
                    //                        {
                    //                            start_index = atoi( token.c_str() );
                    //                        }
                    //                        s.erase( 0, pos + delimiter.length() );
                    //                        counter += 1;
                    //                    }
                    int start_index = Integer.parseInt( tokens[1] );
                    int end_index = Integer.parseInt( tokens[2] );//atoi( s.c_str() );
                    //std::cout + "START INDEX : " + start_index + \n;
                    //std::cout + "END INDEX : " + end_index + \n;
                    phenotype.cycle.data.setTransitionRate( start_index, end_index,
                            phenotype.intracellular.getParameterValue( elm.getValue() ) );
                }
                else
                {
                }
            }

            // volume params
            else if( first.startsWith( "v" ) )
            {
                if( first.equals( "vtsc" ) )
                {
                    phenotype.volume.target_solid_cytoplasmic = phenotype.intracellular.getParameterValue( elm.getValue() );
                }
                else if( first.equals( "vtsn" ) )
                {
                    phenotype.volume.target_solid_nuclear = phenotype.intracellular.getParameterValue( elm.getValue() );
                }
                else if( first.equals( "vff" ) )
                {
                    phenotype.volume.target_fluid_fraction = phenotype.intracellular.getParameterValue( elm.getValue() );
                }
                else
                {
                }
            }
            else
            {
            }

        }
        //std::cout + \n;
        return 0;
    }

    public int validate_PhysiCell_tokens(Microenvironment m, Phenotype phenotype)
    {
        for( Entry<String, String> elm : phenotype_species.entrySet() )
        {
            //std::cout + "PhysiCell_token_validation" + \n;
            //std::cout + elm.first + " : " + elm.second + \n;

            // motility params
            //            String elm.fi
            String first = elm.getKey();
            if( first.startsWith( "m" ) )
            {
                if( first.equals( "mms" ) )
                {
                }
                else if( first.equals( "mpt" ) )
                {
                }
                else if( first.equals( "mmb" ) )
                {
                }
                else
                {
                    //                    System.out.println(" \n;
                    System.out.println( "ERROR: There is no specified token parameters in the name of \"" + first
                            + "\" at motility parameters. Please take a look token specifications." );
                    //                    System.out.println(" \n;
                    //                    System.out.println(" \n;
                    //                    exit (-1);
                    return -1;
                }
            }
            // death params
            else if( first.startsWith( "d" ) )
            {
                if( first.equals( "da" ) )
                {
                }
                else if( first.equals( "dn" ) )
                {
                }
                else
                {
                    //                    System.out.println(" \n;
                    System.out.println( "ERROR: There is no specified token parameters in the name of \"" + first
                            + "\" at death parameters. Please take a look token specifications." );
                    //                    System.out.println(" \n;
                    //                    System.out.println(" \n;
                    //                    exit (-1);
                    return -1;
                }
            }
            // secretion params
            else if( first.startsWith( "s" ) )
            {
                // parsing attribute and getting substrate name
                //                String s = first;
                //                String delimiter = "_";
                //                int pos = 0;
                //                String token;
                //                while( ( pos = s.find( delimiter ) ) != String::npos )
                //                {
                //                    token = s.substr( 0, pos );
                //                    s.erase( 0, pos + delimiter.length() );
                //                }
                String s = first.split( " " )[1];
                int sub_index = m.findDensityIndex( s );
                //std::cout + "SUBSTRATE_INDEX = : " + sub_index + \n;
                if( sub_index < 0 )
                {
                    //                    System.out.println(" \n;
                    System.out.println( "ERROR: There is no substrate named in the name of \"" + s
                            + "\" at microenvironment. Please take a look token specifications." );
                    //                    System.out.println(" \n;
                    //                    System.out.println(" \n;
                    //                    exit (-1);
                    return -1;
                }

                if( first.substring( 0, 3 ).equals( "sur" ) )
                {
                }
                else if( first.substring( 0, 3 ).equals( "ssr" ) )
                {
                }
                else if( first.substring( 0, 3 ).equals( "ssd" ) )
                {
                }
                else if( first.substring( 0, 3 ).equals( "ser" ) )
                {
                }
                else
                {
                    //                    System.out.println(" \n;
                    System.out.println( "ERROR: There is no specified token parameters in the name of \"" + first
                            + "\" at secretion parameters. Please take a look token specifications." );
                    //                    System.out.println(" \n;
                    //                    System.out.println(" \n;
                    //                    exit (-1);
                    return -1;
                }
            }
            else if( first.startsWith( "c" ) )
            {
                if( first.substring( 0, 3 ).equals( "ctr" ) )
                {
                    // getting num of phases
                    int num_of_phases = ( ( phenotype.cycle ) ).phases.size();
                    //std::cout + num_of_phases + \n;

                    // getting start and end indices
                    //                    String s = first;
                    //                    String delimiter = "_";
                    //                    int pos = 0;
                    //                    String token;
                    //                    int counter = 0;
                    //                    int start_index;
                    //                    while( ( pos = s.find( delimiter ) ) != String::npos )
                    //                    {
                    //                        token = s.substr( 0, pos );
                    //                        if( counter == 1 )
                    //                        {
                    //                            start_index = atoi( token.c_str() );
                    //                        }
                    //                        s.erase( 0, pos + delimiter.length() );
                    //                        counter += 1;
                    //                    }
                    String[] tokens = first.split( "_" );
                    int start_index = Integer.parseInt( tokens[1] );
                    int end_index = Integer.parseInt( tokens[2] );

                    // validating the indices
                    if( start_index > num_of_phases - 1 )
                    {
                        //                        System.out.println(" \n;
                        System.out.println( "ERROR: Given transition start index is beyond cycle indices. Please double check it." );
                        //                        System.out.println(" \n;
                        //                        System.out.println(" \n;
                        //                        exit (-1);
                        return -1;
                    }
                    if( end_index > num_of_phases - 1 )
                    {
                        //                        System.out.println(" \n;
                        System.out.println( "ERROR: Given transition end index is beyond cycle indices. Please double check it." );
                        //                        System.out.println(" \n;
                        //                        System.out.println(" \n;
                        //                        exit (-1);
                        return -1;
                    }
                }
                else
                {
                    //                    System.out.println(" \n;
                    System.out.println( "ERROR: There is no specified token parameters in the name of \"" + first
                            + "\" at cycle parameters. Please take a look token specifications." );
                    //                    System.out.println(" \n;
                    //                    System.out.println(" \n;
                    //                    exit (-1);
                    return -1;
                }
            }

            else if( first.startsWith( "v" ) )
            {
                if( first.equals( "vtsc" ) )
                {
                }
                else if( first.equals( "vtsn" ) )
                {
                }
                else if( first.equals( "vff" ) )
                {
                }
                else
                {
                    //                    System.out.println(" \n;
                    System.out.println( "ERROR: There is no specified token parameters in the name of \"" + first
                            + "\" at volume parameters. Please take a look token specifications." );
                    //                    System.out.println(" \n;
                    //                    System.out.println(" \n;
                    //                    exit (-1);
                    return -1;
                }
            }
            else
            {
                //                System.out.println(" \n;
                System.out.println( "ERROR: There is no specified token parameters in the name of \"" + first
                        + "\" at phenotypic parameters. Please take a look token specifications." );
                //                System.out.println(" \n;
                //                System.out.println(" \n;
                //                exit (-1);
                return -1;
            }

        }
        System.out.println( "---- Specified PhysiCell tokens at config file are validated. ----- " );

        return 0;
    }
}