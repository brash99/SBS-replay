// This script is used to add x,y,z offsets to each and every FPP GEM module DB position enetry (sbs.gem<FTorFPP>.position).
// How to use this script:
// 1 - Copy paste the GEM module ang and pos entries from the DB and create a text file.
// 2 - Run this script: '$root 'GEM_addDPposOffset(xOffset, yOffset, zOffset, "<in_filename", "<out_filename>")''
// 3 - Copy and paste the new values in the output file into the DB.

int GEM_addDBposOffset( const double xOffset, const double yOffset, const double zOffset, const std::string in_filename, const std::string out_filename = "outDBGEMPosWithOffsets.txt" )
{
    std::ifstream infile(in_filename);

    if ( !infile.is_open() ) 
    {
        std::cerr << "Eroor: Can't open the input file!" << std::endl;
        
        return -1;
    }

    std::ofstream outfile(out_filename);

    // Read the input file line by line.
    TString currentline;

    while ( currentline.ReadLine(infile) )
    {
        // Ignore comment lines or empty lines
        if ( !currentline.BeginsWith("#") )
        {
            if( currentline.Contains("sbs.gemFPP.m") && currentline.Contains(".position") )
            {
                TPRegexp re("^sbs\\.gemFPP\\.m[0-9]+\\.position");
                if (re.Match(currentline)) // It matches the pattern
                {
                    TObjArray *tokens = currentline.Tokenize(" ");

				    int ntokens = tokens->GetEntries();

                    if ( ntokens >= 5 )
                    {
                        TString xPosStr = ( (TObjString*)(*tokens)[2] )->GetString();
                        TString yPosStr = ( (TObjString*)(*tokens)[3] )->GetString();
                        TString zPosStr = ( (TObjString*)(*tokens)[4] )->GetString();

                        double xPosWithOffset = std::stod(xPosStr.Data()) + xOffset;
                        double yPosWithOffset = std::stod(yPosStr.Data()) + yOffset;
                        double zPosWithOffset = std::stod(zPosStr.Data()) + zOffset;

                        std::setprecision(12);

                        TString newLineWithOffsets = TString::Format("%s = %.10f %.10f %.10f", ((TObjString*)(*tokens)[0])->GetString().Data(), xPosWithOffset, yPosWithOffset, zPosWithOffset);
                        
                        outfile << newLineWithOffsets;
                    }              

                    outfile << std::endl;
                }
            }
            else if ( currentline.Contains("sbs.gemFPP.m") && currentline.Contains(".angle") )
            {
                TPRegexp re("^sbs\\.gemFPP\\.m[0-9]+\\.angle");
                if (re.Match(currentline)) 
                {
                    outfile << currentline;
                    outfile << std::endl;
                }
            }
        }
    }

    infile.close();

    outfile.close();

    std::cout << "New file: " << out_filename << " created with adding the given offsets to GEM module position entries." << std::endl;

    return 0;

}