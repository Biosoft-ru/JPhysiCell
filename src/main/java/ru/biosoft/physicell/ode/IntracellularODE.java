package ru.biosoft.physicell.ode;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.CellDefinition;
import ru.biosoft.physicell.core.Intracellular;
import ru.biosoft.physicell.core.Model;
import ru.biosoft.physicell.core.Phase;
import ru.biosoft.physicell.core.PhaseLink;
import ru.biosoft.physicell.core.Phenotype;

public abstract class IntracellularODE extends Intracellular
{
    protected Map<String, String> phenotypeSpecies = new HashMap<>();

    public IntracellularODE(Model model, CellDefinition cd)
    {
        super( model, cd );
    }

    public Map<String, String> getPhenotypeMapping()
    {
        return phenotypeSpecies;
    }

    public static Map<String, String> createDictionary(String[] densities, CellDefinition cd)
    {
        Map<String, String> dictionary = new TreeMap<>();
        dictionary.put( "Migration speed", "mms" );
        dictionary.put( "Migration persistence time", "mpt" );
        dictionary.put( "Migration bias", "mmb" );

        dictionary.put( "Apoptosis rate", "da" );
        dictionary.put( "Necrosis rate", "dn" );

        for( String density : densities )
        {
            dictionary.put( density + " density", density );
            dictionary.put( density + " uptake rate", density + "_sur" );
            dictionary.put( density + " secretion rate", density + "_ssr" );
            dictionary.put( density + " saturation density", density + "_ssd" );
            dictionary.put( density + " export rate", density + "_ser" );
        }

        dictionary.put( "Target Solid Cytoplasmic Volume", "vtsc" );
        dictionary.put( "Target Solid Nuclear Volume", "vtsn" );
        dictionary.put( "Target Fluid Fraction Volume", "vff" );

        for( List<PhaseLink> links : cd.phenotype.cycle.phaseLinks )
        {
            for( PhaseLink link : links )
            {
                Phase start = link.getStartPhase();
                Phase end = link.getEndPhase();
                dictionary.put( start.name + " -> " + end.name, "ctr_" + start.index + "_" + end.index );
            }
        }
        return dictionary;
    }

    @Override
    public void updatePhenotypeParameters(Microenvironment m, Phenotype phenotype) throws Exception
    {
        Set<String> substrateSet = Stream.of( m.densityNames ).collect( Collectors.toSet() );

        for( String phenotypeName : getOutputs() )
        {
            String variableName = phenotypeSpecies.get( phenotypeName );
            // motility params
            if( phenotypeName.startsWith( "m" ) )
            {
                if( phenotypeName.equals( "mms" ) )
                {
                    phenotype.motility.migrationSpeed = getParameterValue( variableName );
                }
                else if( phenotypeName.equals( "mpt" ) )
                {
                    phenotype.motility.persistenceTime = getParameterValue( variableName );
                }
                else if( phenotypeName.equals( "mmb" ) )
                {
                    phenotype.motility.migrationBias = getParameterValue( variableName );
                }
            }
            // death params
            else if( phenotypeName.startsWith( "d" ) )
            {
                if( phenotypeName.equals( "da" ) )
                {
                    phenotype.death.rates.set( 0, getParameterValue( variableName ) );
                }
                else if( phenotypeName.equals( "dn" ) )
                {
                    phenotype.death.rates.set( 1, getParameterValue( variableName ) );
                }
            }
            // secretion params
            else if( phenotypeName.startsWith( "s" ) )
            {
                String[] tokens = phenotypeName.split( "_" );
                int sub_index = m.findDensityIndex( tokens[1] );

                //uptake rate
                if( phenotypeName.substring( 0, 3 ).equals( "sur" ) )
                {
                    phenotype.secretion.uptakeRates[sub_index] = getParameterValue( variableName );
                }
                //secretion rate
                else if( phenotypeName.substring( 0, 3 ).equals( "ssr" ) )
                {
                    phenotype.secretion.secretionRates[sub_index] = getParameterValue( variableName );
                }
                //secretion density
                else if( phenotypeName.substring( 0, 3 ).equals( "ssd" ) )
                {
                    phenotype.secretion.saturationDensities[sub_index] = getParameterValue( variableName );
                }
                //net export rate
                else if( phenotypeName.substring( 0, 3 ).equals( "ser" ) )
                {
                    phenotype.secretion.netExportRates[sub_index] = getParameterValue( variableName );
                }
            }
            // cycle params
            else if( phenotypeName.startsWith( "c" ) )
            {
                if( phenotypeName.substring( 0, 3 ).equals( "ctr" ) )
                {
                    String delimiter = "_";
                    String[] tokens = phenotypeName.split( delimiter );
                    int start_index = Integer.parseInt( tokens[1] );
                    int end_index = Integer.parseInt( tokens[2] );
                    phenotype.cycle.data.setTransitionRate( start_index, end_index, getParameterValue( variableName ) );
                }
            }
            // volume params
            else if( phenotypeName.startsWith( "v" ) )
            {
                if( phenotypeName.equals( "vtsc" ) )
                {
                    phenotype.volume.target_solid_cytoplasmic = getParameterValue( variableName );
                }
                else if( phenotypeName.equals( "vtsn" ) )
                {
                    phenotype.volume.target_solid_nuclear = getParameterValue( variableName );
                }
                else if( phenotypeName.equals( "vff" ) )
                {
                    phenotype.volume.target_fluid_fraction = getParameterValue( variableName );
                }
            }
            else if( phenotypeName.startsWith( "intracellular" ) )
            {
                String delimiter = "_";
                String[] tokens = phenotypeName.split( delimiter );
                String substrate = tokens[1];
                double value = getParameterValue( variableName );
                int index = m.getSubstrateIndex( substrate );
                phenotype.molecular.internSubstrates[index] = value * phenotype.volume.total;
            }
            else if( substrateSet.contains( phenotypeName ) )
            {
                double value = getParameterValue( variableName );
                int index = m.getSubstrateIndex( phenotypeName );
                phenotype.molecular.internSubstrates[index] = value * phenotype.volume.total;
            }
        }
    }

    @Override
    public void updateIntracellularParameters(Microenvironment m, Phenotype phenotype) throws Exception
    {
        Set<String> substrateSet = Stream.of( m.densityNames ).collect( Collectors.toSet() );

        for( String phenotypeName : getInputs() )
        {
            String variableName = phenotypeSpecies.get( phenotypeName );
            double value = Double.NEGATIVE_INFINITY;
            if( phenotypeName.startsWith( "m" ) ) // motility params
            {
                if( phenotypeName.equals( "mms" ) )
                {
                    value = phenotype.motility.migrationSpeed;
                }
                else if( phenotypeName.equals( "mpt" ) )
                {
                    value = phenotype.motility.persistenceTime = getParameterValue( variableName );
                }
                else if( phenotypeName.equals( "mmb" ) )
                {
                    value = phenotype.motility.migrationBias = getParameterValue( variableName );
                }
            }
            else if( phenotypeName.startsWith( "d" ) ) // death params
            {
                if( phenotypeName.equals( "da" ) )
                {
                    value = phenotype.death.rates.set( 0, getParameterValue( variableName ) );
                }
                else if( phenotypeName.equals( "dn" ) )
                {
                    value = phenotype.death.rates.set( 1, getParameterValue( variableName ) );
                }
            }
            else if( phenotypeName.startsWith( "s" ) ) // secretion params
            {
                String[] tokens = phenotypeName.split( "_" );
                int sub_index = m.findDensityIndex( tokens[1] );
                if( phenotypeName.substring( 0, 3 ).equals( "sur" ) )//uptake rate
                {
                    value = phenotype.secretion.uptakeRates[sub_index] = getParameterValue( variableName );
                }
                else if( phenotypeName.substring( 0, 3 ).equals( "ssr" ) ) //secretion rate
                {
                    value = phenotype.secretion.secretionRates[sub_index] = getParameterValue( variableName );
                }
                else if( phenotypeName.substring( 0, 3 ).equals( "ssd" ) ) //secretion density
                {
                    value = phenotype.secretion.saturationDensities[sub_index] = getParameterValue( variableName );
                }
                else if( phenotypeName.substring( 0, 3 ).equals( "ser" ) ) //net export rate
                {
                    value = phenotype.secretion.netExportRates[sub_index] = getParameterValue( variableName );
                }
            }
            else if( phenotypeName.startsWith( "c" ) ) // cycle params
            {
                if( phenotypeName.substring( 0, 3 ).equals( "ctr" ) )
                {
                    String delimiter = "_";
                    String[] tokens = phenotypeName.split( delimiter );
                    int start_index = Integer.parseInt( tokens[1] );
                    int end_index = Integer.parseInt( tokens[2] );
                    value = phenotype.cycle.data.getTransitionRate( start_index, end_index );
                }
            }

            else if( phenotypeName.startsWith( "v" ) ) // volume params
            {
                if( phenotypeName.equals( "vtsc" ) )
                {
                    value = phenotype.volume.target_solid_cytoplasmic = getParameterValue( variableName );
                }
                else if( phenotypeName.equals( "vtsn" ) )
                {
                    value = phenotype.volume.target_solid_nuclear = getParameterValue( variableName );
                }
                else if( phenotypeName.equals( "vff" ) )
                {
                    value = phenotype.volume.target_fluid_fraction = getParameterValue( variableName );
                }
            }
            else if( phenotypeName.startsWith( "intracellular" ) )
            {
                String delimiter = "_";
                String[] tokens = phenotypeName.split( delimiter );
                String substrate = tokens[1];
                //                double value = getParameterValue( variableName );
                int index = m.getSubstrateIndex( substrate );
                value = phenotype.molecular.internSubstrates[index] / phenotype.volume.total;
                //                phenotype.molecular.internSubstrates[index] = value * phenotype.volume.total;
            }
            else if( substrateSet.contains( phenotypeName ) )
            {
                //                double value = getParameterValue( variableName );
                int index = m.getSubstrateIndex( phenotypeName );
                value = phenotype.molecular.internSubstrates[index] / phenotype.volume.total;
                //                phenotype.molecular.internSubstrates[index] = value * phenotype.volume.total;
            }
            this.setParameterValue( variableName, value );
        }
    }


    public int validate_PhysiCell_tokens(Microenvironment m, Phenotype phenotype)
    {
        for( Entry<String, String> elm : phenotypeSpecies.entrySet() )
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