#include"mixer.h"

Mixer::Mixer(int num_signals, int signal_duration)
{
    mixing_mat = std::make_unique<Mat>(Mat::Zero(num_signals, signal_duration));
    setNumSignalsAndSamples(num_signals,signal_duration);
}

Mixer::Mixer(int num_signals, int signal_duration, Eigen::MatrixXd A)
{
    mixing_mat = std::make_unique<Mat>(Mat::Zero(num_signals, signal_duration));
    setNumSignalsAndSamples(num_signals,signal_duration);
    setMixingMatrix(A);
}

void Mixer::setNumSignalsAndSamples(int new_num_sigs, int new_num_samps)
{
    if(new_num_sigs < 1)
        throw std::invalid_argument("Invalid number of signals, must be greater than 0");
    if(new_num_samps < 1)
        throw std::invalid_argument("Invalid sample size, must be greater than 0");

    n_sig = new_num_sigs;
    n_samp = new_num_samps;
        
    raw_sigs = std::make_shared<Mat>(Mat::Zero(n_sig, n_samp));
    mixing_mat = std::make_unique<Mat>(Mat::Zero(mixing_mat->rows(), n_sig));
    mixed_sigs = std::make_shared<Mat>(Mat::Zero(mixing_mat->rows(), n_samp));
}

void Mixer::setNumSignals(int new_num_sigs)
{

    if(new_num_sigs < 1)
        throw std::invalid_argument("Invalid number of signals, must be greater than 0");

    n_sig = new_num_sigs;
    raw_sigs = std::make_shared<Mat>(Mat::Zero(n_sig, n_samp));
    mixing_mat = std::make_unique<Mat>(Mat::Zero(mixing_mat->rows(), n_sig));
}

void Mixer::setNumSamples(int new_num_samps)
{
    if(new_num_samps < 1)
        throw std::invalid_argument("Invalid sample size, must be greater than 0");

    n_samp = new_num_samps;
    raw_sigs = std::make_shared<Mat>(Mat::Zero(n_sig, n_samp));
    mixed_sigs = std::make_shared<Mat>(Mat::Zero(mixing_mat->rows(), n_samp));
}

void Mixer::setMixingMatrix(Mat A)
{
    if(A.cols() != raw_sigs->rows())
        throw std::invalid_argument("Mixing matrix A's number of columns must match the number of signals to be mixed");

    mixing_mat = std::make_unique<Mat>(A);

    if(mixed_sigs->rows() != A.rows())
        mixed_sigs = std::make_shared<Mat>(Mat::Zero(A.rows(),n_samp));
}


void Mixer::setSourceFromFile(std::string& file_name)
/*Read in raw sources as a matrix from file_name
Format:
number_of_rows number_of_columns T(optional if you want matrix transposed)
x x x ....
.
.
.
where x are doubles
      number_of_rows number_of_columns are ints
*/
{
    std::ifstream mat_stream (file_name);
    if (!mat_stream)
        throw std::invalid_argument("File with name \"" +file_name+ "\" not found");

    std::string line, elem, transposed{"None"};
    bool is_transposed{false};
    int rows{0},cols{0};

    std::getline (mat_stream, line);
    std::stringstream ss(line);

    ss >> rows >> cols;
    if(ss.fail())
        throw std::invalid_argument("Invalid matrix row params: /n"+line);
    if(rows < 1)
        throw std::invalid_argument("Invalid number of rows, must be greater than 0");
    if(cols < 1)
        throw std::invalid_argument("Invalid number of columns, must be greater than 0");

    if(ss>>transposed){
        if(transposed == "T")
            is_transposed = true;
        else
          throw std::invalid_argument("Invalid option. Use T to indicate the matrix must be transposed");  
    }

    Mat SS = Mat::Zero(rows,cols);

    int i{0};
    while (std::getline (mat_stream, line)) {

        std::stringstream ss(line);
        int j{0};
        if(i == rows)
                throw std::invalid_argument("Too many rows:");

        while (std::getline(ss, elem, ','))
        {
            if(j == cols)
                throw std::invalid_argument("Too many elements on row:" +std::to_string(i));

            if(ss.fail())
                throw std::invalid_argument("Invalid matrix row, line: " + std::to_string(i+1)+"/n"+line);
            
            if(!elem.empty())
                SS(i,j++) = stod(elem);
            else
                throw std::invalid_argument("Error on line: " +std::to_string(i+1)+"/n"+line);
        }
        if(j != cols)
            throw std::invalid_argument("Too few elements on row:" +std::to_string(i));
        i++;
    }

    if(i != rows)
        throw std::invalid_argument("Too few rows:");

    mat_stream.close();

    if(!is_transposed){
        setNumSignalsAndSamples(rows, cols);
        raw_sigs = std::make_shared<Mat>(SS);
    }
    else{
        setNumSignalsAndSamples(cols, rows);
        raw_sigs = std::make_shared<Mat>(SS.transpose());
    }
}

void Mixer::setMixingFromFile(std::string& file_name)
/*Read in mixing matrix from file_name
Format:
number_of_rows number_of_columns T(optional if you want matrix transposed)
x x x ....
.
.
.
where x are doubles
      number_of_rows number_of_columns are ints
*/
{
   std::ifstream mat_stream (file_name);
    if (!mat_stream)
        throw std::invalid_argument("File with name \"" +file_name+ "\" not found");

    std::string line, elem, transposed{"None"};
    bool is_transposed{false};
    int rows{0},cols{0};

    std::getline (mat_stream, line);
    std::stringstream ss(line);

    ss >> rows >> cols;
    if(ss.fail())
        throw std::invalid_argument("Invalid matrix row params: /n"+line);

    if(rows < 1)
        throw std::invalid_argument("Invalid number of rows, must be greater than 0");
    if(cols < 1)
        throw std::invalid_argument("Invalid number of columns, must be greater than 0");

    if(ss>>transposed){
        if(transposed == "T")
            is_transposed = true;
        else
          throw std::invalid_argument("Invalid option. Use T to indicate the matrix must be transposed");  
    }
    
    Mat AA = Mat::Zero(rows,cols);

    int i{0};
    while (std::getline (mat_stream, line)) {

        std::stringstream ss(line);
        int j{0};
        if(i == rows)
                throw std::invalid_argument("Too many rows:");

        while (std::getline(ss, elem, ','))
        {
            if(j == cols)
                throw std::invalid_argument("Too many elements on row:" +std::to_string(i));

            if(ss.fail())
                throw std::invalid_argument("Invalid matrix row, line: " + std::to_string(i+1)+"/n"+line);
            
            if(!elem.empty())
                AA(i,j++) = stod(elem);
            else
                throw std::invalid_argument("Error on line: " +std::to_string(i+1)+"/n"+line);
        }
        if(j != cols)
            throw std::invalid_argument("Too few elements on row:" +std::to_string(i));
        i++;
    }

    if(i != rows)
        throw std::invalid_argument("Too few rows:");
    
    mat_stream.close();

    if(!is_transposed)
        setMixingMatrix(AA);
    else
        setMixingMatrix(AA.transpose());

}

void Mixer::genSignals(std::string& gen_file )
/* Textfile must be of form

number_of_signals number_of_samples
signal_amplitude_1 signal_period_1 signal_type_1
signal_amplitude_2 signal_period_2 signal_type_2

where:
    signal_type is a string
    number_of_signals and number_of_samples are integers > 0
    signal_amplitude and signal_period are double
*/
{
    std::ifstream instr_file (gen_file);
    if(!instr_file)
        throw std::invalid_argument("File with name \"" +gen_file+ "\" not found");

    std::string line,opt;
    int m,n,i{0};
    double amp, tau;

    std::getline (instr_file, line);
    std::stringstream ss(line);
    ss>>m>>n;

    if(ss.fail())
        throw std::invalid_argument("Invalid instruction in text file, line: " + std::to_string(i+2));

    setNumSignalsAndSamples(m,n);
    SigGen::WaveGen waves((int)n_samp);

    while (std::getline (instr_file, line)) {

        std::stringstream ss(line);
        ss >> amp >> tau >> opt;

        if(ss.fail())
            throw std::invalid_argument("Invalid instruction in text file, line: " + std::to_string(i+2));

        raw_sigs->row(i++) = waves.genWave(amp, tau, opt);

        if(i > n_sig)
            throw std::invalid_argument("Mismatch between number of signals (" +std::to_string(n_sig)+ 
                ") and amount of signal instructions in file (" +std::to_string(i)+")");
    }

    instr_file.close();

}

void Mixer::genSignals(int seed, int max_amp, int max_period)
{
    srand(seed);
    SigGen::WaveGen waves((int)n_samp);

    for(int i=0; i<n_sig; i++)
    {
        raw_sigs->row(i) = waves.genRandomWave(max_amp,max_period); //default max_amp=2 max_period=100
    }

}

void Mixer::mixSignals(bool noisy,int seed)
{
    if(noisy)
    {
        std::mt19937 r_gen (seed);
        std::normal_distribution<double> r_dis(0.0, 0.1);
        auto norm = [&] () {return r_dis(r_gen);};

        Mat v = Mat::NullaryExpr(n_sig,n_samp, norm );

        *raw_sigs += v;
        raw_sigs->rowwise().normalize();
    }
    else
    {
        raw_sigs->rowwise().normalize();
    }

    *mixed_sigs = *mixing_mat * *raw_sigs;
}
