#include "src/XamReader.h"
#include "src/XamTag.h"
#include <vector>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
    // gather file names
    vector<string> filenames;
    for (int i = 1; i < argc; i++)
    {
        filenames.push_back(argv[i]);
    }


    cout << "filenum\tchrom\tstart\tend\tname\tisize\tis_discordant\tread_group\tedit_dist\talign_score\tbases\n";
    XamReader *xr = new XamReader(filenames);
    while (xr->next() != false)
    {
        cout << xr->rec->filenum
             << "\t"
             << xr->rec->chrom()
             << "\t"
             << xr->rec->alignment_start()
             << "\t"
             << xr->rec->alignment_stop()
             << "\t"
             << xr->rec->name()
             << "\t"
             << xr->rec->insert_size()
             << "\t"
             << !(xr->rec->properly_paired())
             << "\t"
             << xr->rec->read_group()
             << "\t"
             << xr->rec->edit_dist()
             << "\t"
             << xr->rec->alignment_score()
             << "\t"
             << xr->rec->bases()
             << endl;
        delete xr->rec;
    }
}