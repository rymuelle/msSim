#include <iostream>
#include <TROOT.h>
#include "math.h"
#include <random>
#include <string>

// g++ -o msSim -c msSim.C
/*
def createScatteringRanges(low,high):

    numScatter =  random.randint(low,high)
    begin, end = 0, 1
    for j in range(numScatter):
        msScatterDistance.append(random.uniform(begin,end))
        #begin = msScatterDistance[j]
    #print msScatterDistance
    #print random.gauss(0, 0.1)
    msScatterDistance.sort()
*/
std::vector<double> createScatteringRanges(std::exponential_distribution<> rng, std::mt19937 & rnd_gen , double pathLength){

    std::vector<double> scatter_distances;
    double scatter_distance;
    double total_distance = 0;
    while(total_distance < pathLength*1.5){
        scatter_distance = (double)rng (rnd_gen) * 4;
        
        scatter_distances.push_back(scatter_distance);

        total_distance = total_distance + scatter_distance;
    
    }

    return scatter_distances;
}

std::vector<double> createScatteringAngles(std::vector<double> scatter_vector, double scatter_angle, std::mt19937 & rnd_gen, std::normal_distribution<double> scatter_dist){

    std::vector<double> scattering_angles_vector;

    for(int i =0; i< scattering_angles_vector.size(); i++){

         std::cout << scatter_dist(rnd_gen) << std::endl;

    }

    return scattering_angles_vector;
}


TH1F runSim(double x_pos_detector, double y_target_pos, int num, double scatter_angle, std::mt19937 & rnd_gen, std::normal_distribution<double> scatter_angle_dist, std::exponential_distribution<> rng){

   
    TH1F TH1F_x_residual = TH1F("TH1F_x_residual", "x residual at ideal position; x; counts", 100, -1,1);

    for(int i=0;i<num;i++){
/*
        scatter_vector = createScatteringRanges(rng,rnd_gen, pathLength);
        createScatteringAngles(scatter_vector, scatter_angle, rnd_gen, scatter_dist);
*/
        

        std::vector<double> scatter_distances;
        double scatter_distance;
        double total_distance = 0;
        std::vector<double> scattering_angles_vector;
        double x_pos = 0;
        double y_pos = 0;

        std::vector<double> x_pos_vec;
        std::vector<double> y_pos_vec;


        double x_distance;
        double y_distance;
        double scatter_angle_l;

        double angle = atan(y_target_pos/x_pos_detector);
        while(true){
            scatter_distance = (double) rng (rnd_gen) * 4;

            scatter_angle_l = scatter_angle_dist(rnd_gen);

            //std::cout << scatter_angle_l << " " << scatter_distance << " " << x_pos << " "  << y_pos  << std::endl;
    
            total_distance = total_distance + scatter_distance;

            x_pos = scatter_distance*cos(angle) + x_pos;
            y_pos = scatter_distance*sin(angle) + y_pos;

            angle = angle + scatter_angle_l;

            x_pos_vec.push_back(x_pos);
            y_pos_vec.push_back(y_pos);

            if(x_pos > x_pos_detector){
                x_distance = x_pos_detector - x_pos_vec.at(x_pos_vec.size() -2);

                y_distance = x_distance*y_pos/x_pos + y_pos_vec.at(y_pos_vec.size() -2);

                x_distance = x_distance + x_pos_vec.at(x_pos_vec.size() -2);

                break;
            }

        }
        TH1F_x_residual.Fill(y_distance- y_target_pos);
    
        //std::cout << TH1F_x_residual->GetEntries() << std::endl;

        //std::cout << scatter_angle_l << " " << scatter_distance << " " << x_pos << " "  << y_pos  <<" " << x_distance << " " << y_distance<<  " " << y_distance- y_target_pos << std::endl;

        
        
        if(i*10 % num == 0){
            std::cout <<  (double)i/(double)num << std::endl;
            
        }
    }

    std::cout << "RMS " << TH1F_x_residual.GetRMS() << " mean " << TH1F_x_residual.GetMean() << std::endl;
    TH1F_x_residual.Draw();
    //c1->SaveAs( (std::string("output/TH1F_x_residual_xpos_") + std::to_string(x_pos_detector) +std::string("_ypos_") + std::to_string(y_target_pos)  + std::string(".png" )).c_str() ) ;


    return TH1F_x_residual;

}




int msSim(){
    double x_pos_detector = 1.0;
    double y_target_pos  = 1.0;
    double pathLength = sqrt(x_pos_detector*x_pos_detector + y_target_pos*y_target_pos);

    int num = 100000;
    TCanvas c1;
   // TH1F * TH1F_x_residual = new TH1F("TH1F_x_residual", "x residual at ideal position; x; counts", 100, -1,1);
    TH1F TH1F_x_residual_corrected = TH1F("TH1F_x_residual_corrected", "x residual at ideal position corrected; x; counts", 100, -.5,.5);
    TH2F TH2F_x_residual_dx_residual = TH2F("TH2F_x_residual_dx_residual", "x vs dx/dy residual", 100, -1,1,100, 1,1);
    TH2F TH2F_x_residual_v_ms = TH2F("TH2F_x_residual_v_ms", "x vs dx/dy residual", 100, -1,1,100, 0,100);


    double scatterLengthScale = x_pos_detector/80;
    double scatterLengthScaleLambda = 1/scatterLengthScale;
    double scatter_angle = .0001;

    std::random_device rd; 
    std::mt19937 rnd_gen (rd());

    std::exponential_distribution<> rng (scatterLengthScaleLambda);
    std::normal_distribution<double> scatter_angle_dist(0.0,scatter_angle);

    
    std::vector<double> scatter_vector;

   
    TH1F *Temp = (TH1F*) runSim( x_pos_detector,  y_target_pos,  num,  scatter_angle, rnd_gen, scatter_angle_dist,  rng).Clone();
    Temp->Draw();
    c1.SaveAs( (std::string("output/TH1F_x_residual_xpos_") + std::to_string(x_pos_detector) +std::string("_ypos_") + std::to_string(y_target_pos)  + std::string(".png" )).c_str() ) ;


    //for(int i=0;i<num;i++){
/*
        scatter_vector = createScatteringRanges(rng,rnd_gen, pathLength);
        createScatteringAngles(scatter_vector, scatter_angle, rnd_gen, scatter_dist);
*/
/*
        std::vector<double> scatter_distances;
        double scatter_distance;
        double total_distance = 0;
        std::vector<double> scattering_angles_vector;
        double x_pos = 0;
        double y_pos = 0;

        std::vector<double> x_pos_vec;
        std::vector<double> y_pos_vec;


        double x_distance;
        double y_distance;
        double scatter_angle_l;

        double angle = atan(y_target_pos/x_pos_detector);
        while(true){
            scatter_distance = (double) rng (rnd_gen) * 4;

            scatter_angle_l = scatter_angle_dist(rnd_gen);

            //std::cout << scatter_angle_l << " " << scatter_distance << " " << x_pos << " "  << y_pos  << std::endl;
    
            total_distance = total_distance + scatter_distance;

            x_pos = scatter_distance*cos(angle) + x_pos;
            y_pos = scatter_distance*sin(angle) + y_pos;

            angle = angle + scatter_angle_l;

            x_pos_vec.push_back(x_pos);
            y_pos_vec.push_back(y_pos);

            if(x_pos > x_pos_detector){
                x_distance = x_pos_detector - x_pos_vec.at(x_pos_vec.size() -2);

                y_distance = x_distance*y_pos/x_pos + y_pos_vec.at(y_pos_vec.size() -2);

                x_distance = x_distance + x_pos_vec.at(x_pos_vec.size() -2);

                break;
            }

        }
        TH1F_x_residual->Fill(y_distance- y_target_pos);
    
        //std::cout << TH1F_x_residual->GetEntries() << std::endl;

        std::cout << scatter_angle_l << " " << scatter_distance << " " << x_pos << " "  << y_pos  <<" " << x_distance << " " << y_distance<<  " " << y_distance- y_target_pos << std::endl;

        
        
        if(i*10 % num == 0){
            std::cout <<  (double)i/(double)num << std::endl;
            
        }
    }
   

        std::cout << "RMS " << TH1F_x_residual->GetRMS() << " mean " << TH1F_x_residual->GetMean() << std::endl;
        TH1F_x_residual->Draw();
        
        */
         return 0;



}



