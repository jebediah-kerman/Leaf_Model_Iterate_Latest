
// add to freefem++
// https://doc.freefem.org/documentation/developers.html#adding-a-new-finite-element



 // void init(){
 //     Global.Add("datetime", "(", new OneOperator0s<string>(datetime));
 // }
 // LOADFUNC(init);



#include <iostream>
#include <ctime>

// add datetime for saving files
//https://www.tutorialspoint.com/cplusplus/cpp_date_time.htm

using namespace std;


int main() {
   // current date/time based on current system
   time_t now = time(0);
   
   // convert now to string form
   char* dt = ctime(&now);

   cout << "The local date and time is: " << dt << endl;

   // convert now to tm struct for UTC
   tm *gmtm = gmtime(&now);
   dt = asctime(gmtm);
   cout << "The UTC date and time is:"<< dt << endl;

   // return dt

   cout << GetFolder("AR") << endl;
}




