package ru.biosoft.physicell.ode;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import ru.biosoft.physicell.biofvm.Microenvironment;
import ru.biosoft.physicell.core.Cell;
import ru.biosoft.physicell.core.Intracellular;
import ru.biosoft.physicell.core.Phenotype;
import ru.biosoft.physicell.sample_projects.ode_energy.ToyMetabolicModel;

public class IntracellularEuler extends Intracellular
{

    //    Map<String, Double> parameters;
    //    Map<String, String> substrate_species;
    //    Map<String, String> custom_data_species;
    protected Map<String, String> phenotype_species = new HashMap<>();
    //    Map<String, Integer> speciesIndex;
    double next_librr_run = 0;

    private ODEModel model;
    private Euler euler;


    public IntracellularEuler() throws Exception
    {
        this.model = new ToyMetabolicModel();
        euler = new Euler( model );
    }

    public IntracellularEuler(ODEModel model) throws Exception
    {
        this.model = model;
        euler = new Euler( model );
        //        intracellular_type = "sbml";
        //        std::cout + "====== " + __FUNCTION __ + "() intracellular_type=" + intracellular_type + \n;
        //        std::cout + "====== " + __FUNCTION__ + "() sbml_filename = " +  sbml_filename + \n;
        // initial_values.clear();
        // mutations.clear();
        //        parameters.clear();
    }

    // constructor using XML node
    //    RoadRunnerIntracellular::RoadRunnerIntracellular(pugi::xml_node& node)
    //    {
    // std::cout + "======rwh " + __FUNCTION__ + ": node.name() =" + node.name() + \n;
    //        intracellular_type = "roadrunner";
    //        initialize_intracellular_from_pugixml(node);
    // std::cout + "======rwh " + __FUNCTION__ + "(node) intracellular_type=" + intracellular_type + \n;
    // std::cout + "======rwh " + __FUNCTION__ + "(node) sbml_filename = " +  sbml_filename + \n;
    // std::cout + "======rwh " + __FUNCTION__ + "(node) this=" +  this + \n;
    // std::cout + "======rwh " + __FUNCTION__ + "(node) this.sbml_filename=" +  this.sbml_filename + \n;
    //    }

    // Intracellular* RoadRunnerIntracellular::clone() // -. 'Intracellular' does not name a type
    // {
    //  return static_cast<Intracellular*>(new RoadRunnerIntracellular(this));
    // }

    // rwh: review this
    //    RoadRunnerIntracellular::RoadRunnerIntracellular(RoadRunnerIntracellular* copy) 
    //    {
    //        intracellular_type = copy.intracellular_type;
    //        sbml_filename = copy.sbml_filename;
    // cfg_filename = copy.cfg_filename;
    // time_step = copy.time_step;
    // discrete_time = copy.discrete_time;
    // time_tick = copy.time_tick;
    // scaling = copy.scaling;
    // initial_values = copy.initial_values;
    // mutations = copy.mutations;
    //        parameters = copy.parameters;
    //        
    //    }

    // Parse the <intracellular> info in the .xml for (possibly) each <cell_definition ...>, e.g.
    // <intracellular type="roadrunner">
    //  <sbml_filename>./config/Toy_SBML_Model_2.xml</sbml_filename>
    //  <time_step>1</time_step>
    //         <species substrate="oxygen">Oxy</species>
    //         <species substrate="glucose">Glc</species>
    //         <species custom_data="energy">Energy</species>
    //    void RoadRunnerIntracellular::initialize_intracellular_from_pugixml(pugi::xml_node& node)
    //    {
    //        pugi::xml_node node_sbml = node.child( "sbml_filename" );
    //        if ( node_sbml )
    //        { 
    //            sbml_filename = PhysiCell::xml_get_my_string_value (node_sbml); 
    //            std::cout + "\n------------- "  + __FUNCTION__ + ": sbml_filename = " + sbml_filename + \n;
    //        }
    //        
    //        pugi::xml_node node_species = node.child( "map" );
    //        while( node_species )
    //        {
    //            // ---------  substrate
    //            
    //            String substrate_name = node_species.attribute( "PC_substrate" ).value(); 
    //            if( substrate_name != "" )
    //            {
    //                //std::cout + "-----------" + node_species.attribute( "sbml_species" ).value() + \n; 
    //                String species_name = node_species.attribute( "sbml_species" ).value();
    //                substrate_species[substrate_name] = species_name;
    //                std::cout + "\n------------- "  + __FUNCTION__ + ": species_name= " + species_name + \n;
    //            }
    //            // ---------  custom_data
    //            String custom_data_name = node_species.attribute( "PC_custom_data" ).value(); 
    //            if( custom_data_name != "" )
    //            {
    //                String species_name = node_species.attribute( "sbml_species" ).value();
    //                custom_data_species[custom_data_name] = species_name;
    //                // std::cout + "\n------------- "  + __FUNCTION__ + ": species_name= " + species_name + \n;
    //            }
    //            
    //            
    //            // ---------  phenotype_data
    //            String phenotype_name = node_species.attribute( "PC_phenotype" ).value(); 
    //            
    //            if( phenotype_name != "" )
    //            {
    //                String species_name = node_species.attribute( "sbml_species" ).value();
    //                phenotype_species[phenotype_name] = species_name;
    //                // std::cout + "\n------------- "  + __FUNCTION__ + ": species_name= " + species_name + \n;
    //            }
    //
    //            node_species = node_species.next_sibling( "map" ); 
    //        }
    //        
    //        std::cout + "  ------- substrate_species map:"  + \n;
    //        for(auto elm : substrate_species)
    //        {
    //            std::cout + "      "  + elm.first + " . " + elm.second + \n;
    //        }
    //        std::cout + "  ------- custom_data_species map:"  + \n;
    //        for(auto elm : custom_data_species)
    //        {
    //            std::cout + "      "  + elm.first + " . " + elm.second + \n;
    //        }
    //        std::cout + \n;
    //
    //        std::cout + "  ------- phenotype_species map:"  + \n;
    //        for(auto elm : phenotype_species)
    //        {
    //            std::cout + "      "  + elm.first + " . " + elm.second + \n;
    //        }
    //        std::cout + \n;
    //
    //    }
    //
    //
    //    void RoadRunnerIntracellular::start()
    //    {
    //        // called when a new cell is created; creates the unique 'rrHandle'
    //        rrc::RRVectorPtr vptr;
    //
    //        //std::cout + "\n------------ " + __FUNCTION__ + ": librr_intracellular.cpp: start() called\n";
    //        // this.enabled = true;
    //
    //        //std::cout + "\n------------ " + __FUNCTION__ + ": doing: rrHandle = createRRInstance()\n";
    //
    //        rrHandle = createRRInstance();
    //
    //        //std::cout + "\n------------ " + __FUNCTION__ + ": rrHandle = " + rrHandle + \n;
    //
    //        // if (!rrc::loadSBML (rrHandle, get_cell_definition("lung epithelium").sbml_filename.c_str())) 
    //        //std::cout + "     sbml_filename = " + sbml_filename + \n;
    //
    //        // TODO: don't hard-code name
    //        if ( !rrc::loadSBML(rrHandle, (sbml_filename).c_str() ) )
    //        // std::cout + "     for now, hard-coding sbml_file = ./config/Toy_SBML_Model_1.xml" + \n;
    //        // if (!rrc::loadSBML(rrHandle, "./config/Toy_SBML_Model_1.xml") )
    //        {
    //            std::cerr + "------------.>>>>  Error while loading SBML file  <-------------\n\n";
    //            // return -1;
    //            //  printf ("Error message: %s\n", getLastError());
    //            exit (0);
    //        }
    //
    //        // std::cout + "     rrHandle=" + rrHandle + \n;
    //
    //        int r = rrc::getNumberOfReactions(rrHandle);
    //        int m = rrc::getNumberOfFloatingSpecies(rrHandle);
    //        int b = rrc::getNumberOfBoundarySpecies(rrHandle);
    //        int p = rrc::getNumberOfGlobalParameters(rrHandle);
    //        int c = rrc::getNumberOfCompartments(rrHandle);
    //
    //
    //        //std::cerr + "Number of reactions = " + r + \n;
    //        //std::cerr + "Number of floating species = " + m + \n;  // 4
    //        //std::cerr + "Number of boundary species = " + b + \n;  // 0
    //        //std::cerr + "Number of compartments = " + c + \n;  // 1
    //
    //        //std::cerr + "Floating species names:\n";
    //        //std::cerr + "-----------------------\n";
    //        String species_names_str = stringArrayToString(rrc::getFloatingSpeciesIds(rrHandle));
    //        //std::cerr +  species_names_str +"\n"+ \n;
    //    Stringstream iss(species_names_str);
    //        String species_name;
    //        int idx = 0;
    //        while (iss >> species_name)
    //        {
    //            species_result_column_index[species_name] = idx;
    //            //std::cout + species_name + " . " + idx + \n;
    //            idx++;
    //        }
    //
    //        vptr = rrc::getFloatingSpeciesConcentrations(rrHandle);
    //        //std::cerr + vptr.Count + \n;
    //    /*     for (int kdx=0; kdx<vptr.Count; kdx++)
    //        {
    //            std::cerr + kdx + ") " + vptr.Data[kdx] + \n;
    //        } */
    //        //std::cerr + "----------  end start() -------------\n";
    //        
    //        rrc::freeVector(vptr);
    //        // return 0;
    //    }

    @Override
    public boolean need_update(double curTime)
    {
        return curTime >= this.next_librr_run;
    }

    // solve the intracellular model
    @Override
    public void update() throws Exception
    {
        //        double[] before = model.x_values.clone();
        euler.doStep();
        //        double[] after = model.x_values;
        //        System.out.println( VectorUtil.print( before, 5 ) + " -> " + VectorUtil.print( after, 5 ) );
        //        double start_time = 0.0;
        //        double end_time = 0.01;
        // static double end_time = 10.0;
        // static int num_vals = 1;
        // static int num_vals = 10;
        //        int num_vals = 2;

        // result = rrc::simulateEx (pCell.phenotype.molecular.model_rr, 0, 10, 10);  // start time, end time, and number of points
        //std::cout + __FUNCTION__ + " ----- update(): this=" + this + \n;
        //std::cout + __FUNCTION__ + " ----- update(): rrHandle=" + this.rrHandle + \n;

        // if (this.result != 0)   // apparently not necessary (done in freeRRCData hopefully)
        //    rrc::freeRRCData (this.result);

        //TODO: SIMULATION HERE
        //        this.result = rrc::simulateEx (this.rrHandle, start_time, end_time, num_vals);  // start time, end time, and number of points


        // this.next_librr_run += this.rrHandle.get_time_to_update();
        // std::cout + "----- update(): result=" + result + \n;
        //std::cout + "----- update(): result.CSize=" + this.result.CSize + \n;
        //std::cout + "----- update(): result.RSize=" + this.result.RSize + \n;  // should be = num_vals
        // std::cout + "----- update(): result.ColumnHeaders[0]=" + result.ColumnHeaders[0] + \n;  // = "time"

        // debug - does it generate expected data?
        int index = 0;
        // Print out column headers... typically time and species.
        //        for( int col = 0; col < this.result.CSize; col++ )
        //        {
        // std::cout + result.ColumnHeaders[index++];
        // std::cout + std::left + std::setw(15) + result.ColumnHeaders[index++];
        //std::cout + std::left + this.result.ColumnHeaders[index++];
        // if (col < result.CSize - 1)
        // {
        //  // std::cout + "\t";
        //  std::cout + "  ";
        // }
        //        }
        //std::cout + "\n";

        index = 0;
        // Print out the data
        //        for( int row = 0; row < this.result.RSize; row++ )
        //        {
        //            for( int col = 0; col < this.result.CSize; col++ )
        //            {
        // std::cout + result.Data[index++];
        //std::cout + std::left + std::setw(15) + this.result.Data[index++];
        // if (col < result.CSize -1)
        // {
        //  // std::cout + "\t";
        //  std::cout + "  ";
        // }
        //            }
        // std::cout + "\n";
        //        }
        // int idx = (result.RSize - 1) * result.CSize + 1;
        // std::cout + "Saving last energy value (cell custom var) = " + result.Data[idx] + \n;
        // pCell.custom_data[energy_cell_idx]  = result.Data[idx];

        // return 0;
    }

    @Override
    public double getParameterValue(String param_name) throws Exception
    {
        return model.getCurrentValue( param_name );
        //        rrc::RRVectorPtr vptr;

        //std::cout + "-----------"  + __FUNCTION__ + "-----------" + \n;
        // std::cout + "    substrate_name = " + substrate_name + \n;
        //std::cout + "    param_name = " + param_name + \n;

        // TODO: optimize this eventually
        // std::map<String, int> species_result_column_index;
        // int num_columns = result.CSize;
        // int offset = (num_rows_result_table-1) * result.CSize - 1;
        // int offset = (num_rows_result_table-1) * result.CSize;
        // offset += (num_rows_result_table-1) * result.CSize - 1;

        // int offset = species_result_column_index[name];
        // String species_name = this.substrate_species[substrate_name];
        // std::cout + "    species_name = " + species_name + \n;

        //        TODO: GET RESULTS HERE
        //        vptr = rrc::getFloatingSpeciesConcentrations(this.rrHandle);
        //std::cerr + vptr.Count + \n;
        //        for (int kdx=0; kdx<vptr.Count; kdx++)
        //        {
        //std::cerr + kdx + ") " + vptr.Data[kdx] + \n;
        //        }

        //        int offset = species_result_column_index[param_name];
        //std::cout + "    result offset = "+ offset + \n;
        // double res = this.result.Data[offset];
        //        double res = vptr.Data[offset];
        //std::cout + "    res = " + res + \n;
        //        rrc::freeVector(vptr);
        //        return 0;//res;
    }

    // rwh: might consider doing a multi-[species_name, value] "set" method

    @Override
    public void setParameterValue(String species_name, double value) throws Exception
    {
        model.setCurrentValue( species_name, value );
        //        rrc::RRVectorPtr vptr;
        //
        //        vptr = rrc::getFloatingSpeciesConcentrations(this.rrHandle);
        //        int idx = species_result_column_index[species_name];
        //        vptr.Data[idx] = value;
        //        // rrc::setFloatingSpeciesConcentrations(pCell.phenotype.molecular.model_rr, vptr);
        //        rrc::setFloatingSpeciesConcentrations(this.rrHandle, vptr);
        //        rrc::freeVector(vptr);
        // return 0;
    }


    //    getRoadRunnerModel(Phenotype phenotype) {
    //        return static_cast<RoadRunnerIntracellular*>(phenotype.intracellular);
    //    }

    //    void save_libRR(String path, String index)
    //    {
    //        String state_file_name = path + "/states_" + index + ".dat";
    //    //  String state_file_name = path + "/states_" + index + ".csv";
    //        std::ofstream state_file( state_file_name );
    //        state_file + "---------  dummy output from save_libRR  ---------" + \n;
    //        state_file + "ID,state" + \n;
    //        for( auto cell : *PhysiCell::all_cells )
    //            state_file + cell.ID + "," + cell.phenotype.intracellular.get_state() + \n;
    //        state_file.close();
    //    }

    //    String get_state()
    //    {
    //        return sbml_filename;
    //    }

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
                if( first.substring( 0, 3 ) == "sur" )
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

    @Override
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

                if( first.substring( 0, 3 ) == "sur" )
                {
                }
                else if( first.substring( 0, 3 ) == "ssr" )
                {
                }
                else if( first.substring( 0, 3 ) == "ssd" )
                {
                }
                else if( first.substring( 0, 3 ) == "ser" )
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
                if( first.substring( 0, 3 ) == "ctr" )
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

    //    int validate_SBML_species()
    //    {
    //        //std::cout + "---------VALIDATING_SBML_SPECIES START-------" + \n;
    //        
    //        // reading SBML
    //        rrHandle = createRRInstance();
    //        if ( !rrc::loadSBML(rrHandle, (sbml_filename).c_str() ) )
    //        {
    //            System.out.println( "------------.>>>>  Error while loading SBML file  <-------------\n\n");
    //            return -1;
    //        } 
    //        // getting Species Names
    //        String species_names_str = stringArrayToString(rrc::getFloatingSpeciesIds(rrHandle));
    //    Stringstream iss(species_names_str);
    //        String species_name;
    //        
    //        std::vector<String> all_species {};
    //        
    //        int idx = 0;
    //        while (iss >> species_name)
    //        {
    //            species_result_column_index[species_name] = idx;
    //            all_species.push_back(species_name);
    //            //std::cout + species_name + " . " + idx + \n;
    //            idx++;
    //        }
    //
    //        // Phenotype Species 
    //        for (EntrySet<String, String> elm : phenotype_species.entrySet())
    //        {
    //            boolean exist = 0;
    //           // std::cout + species_name.size() + \n;
    //            for (int i=0; i < all_species.size(); i++)
    //            {
    //                //std::cout + all_species[i] + \n;;
    //                //std::cout + "Comparing " + all_species[i] + " with " + elm.second + \n;
    //                if ( all_species[i] == elm.second )
    //                {
    //                   //std::cout + "And they are the same..... " +\n;
    //                   exist = 1; 
    //                }
    //                idx++;  
    //            }
    //            if (!exist)
    //            {
    ////                System.out.println(" \n;
    //                System.out.println("ERROR: The specified SBML species in the name of \"" + elm.second + "\" at phenotypic species. Please take a look SBML species specifications.");
    ////                System.out.println(" \n;
    ////                System.out.println(" \n;
    ////                exit (-1);
    //                return -1;
    //            }
    //            //std::cout + "existence check : " + elm.second +": " + exist + \n;
    //        }
    //        
    //        // Substrate Species
    //        for (EntrySet<String, String> elm : substrate_species.entrySet())
    //        {
    //            boolean exist = 0;
    //           // std::cout + species_name.size() + \n;
    //            for (int i=0; i < all_species.size(); i++)
    //            {
    //                //std::cout + all_species[i] + \n;;
    //                //std::cout + "Comparing " + all_species[i] + " with " + elm.second + \n;
    //                if ( all_species[i] == elm.second )
    //                {
    //                   //std::cout + "And they are the same..... " +\n;
    //                   exist = 1; 
    //                }
    //                idx++;  
    //            }
    //            if (!exist)
    //            {
    ////                System.out.println(" \n;
    //                System.out.println("ERROR: The specified SBML species in the name of \"" + elm.second + "\" at substrate species. Please take a look SBML species specifications.");
    ////                System.out.println(" \n;
    ////                System.out.println(" \n;
    ////                exit (-1);
    //                return -1;
    //            }
    //            //std::cout + "existence check : " + elm.second +": " + exist + \n;
    //        }    
    //
    //        // custom data species
    //        for (EntrySet<String, String> elm : custom_data_species.entrySet())
    //        {
    //            boolean exist = 0;
    //           // std::cout + species_name.size() + \n;
    //            for (int i=0; i < all_species.size(); i++)
    //            {
    //                //std::cout + all_species[i] + \n;;
    //                //std::cout + "Comparing " + all_species[i] + " with " + elm.second + \n;
    //                if ( all_species[i] == elm.second )
    //                {
    //                   //std::cout + "And they are the same..... " +\n;
    //                   exist = 1; 
    //                }
    //                idx++;  
    //            }
    //            if (!exist)
    //            {
    ////                System.out.println(" \n");
    //                System.out.println("ERROR: The specified SBML species in the name of \"" + elm.second + "\" at substrate species. Please take a look SBML species specifications.");
    ////                System.out.println(" \n");
    ////                System.out.println(" \n");
    ////                exit (-1);
    //                return -1;
    //            }
    //            //std::cout + "existence check : " + elm.second +": " + exist + \n;
    //        }    
    //        
    //        
    //        //std::cout + "---------VALIDATING_SBML_SPECIES END-------" + \n;  
    //        
    //        
    //        System.out.println("---- Specified SBML species at config file are validated. -----");
    //        return 0;
    //    }

    @Override
    public void update(Cell cell, Phenotype phenotype, double dt)
    {
        // TODO Auto-generated method stub

    }

    @Override
    public void inherit(Cell cell)
    {
        // TODO Auto-generated method stub

    }

    @Override
    public void setDT(double dt)
    {
        this.euler.dt = dt;
    }

    @Override
    public void addPhenotypeSpecies(String code, String species)
    {
        this.phenotype_species.put( code, species );
    }

    //    int create_custom_data_for_SBML(Phenotype phenotype)
    //    {
    //std::cout + "Test" + \n;
    //        
    //        return 0; 
    //    }

    @Override
    public Intracellular clone()
    {
        try
        {
            IntracellularEuler result = (IntracellularEuler)super.clone();
            result.phenotype_species = new HashMap<>( phenotype_species );
            result.model = (ODEModel)model.clone();
            result.euler = (Euler)euler.clone( result.model );
            return result;
            
        }
        catch( Exception e )
        {
            throw ( new InternalError( e ) );
        }
    }
}
