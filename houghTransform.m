%**************************************************************************
%* Moba - Hough-Transformation                                            *
%**************************************************************************

%instead of this, call the create_image_with_lines function in the command
%window.
create_image_with_lines('Hausarbeit 1','moba_HA1_Abb3_NZroads.jpg',1);
create_image_with_lines('Hausarbeit 1','moba_HA1_Abb1_Rechteck.bmp',1);
create_image_with_lines('Hausarbeit 1','moba_HA1_Abb2_Scan.bmp',1);


%Main function
%creates an image in the mat-file folder. this image is the read image with 
%the detected lines in it. if there is a border. the border will be cropped. 
%the created image and the hough space will be shown to the user. 
%this feature can be turned off by setting the do_show_image from 1 to 0 in 
%the end of the upper commands.
function create_image_with_lines(path,image_name,do_show_images)

    %reads the image
    picture_name=strcat(strcat(path,'\'),image_name);
    A=imread(picture_name);
    A=im2bw(A);
    %deletes any frame that is around the image
    A = delete_frame (A);
    %gets the accumulator matrix, representation of the image in hough 
    %space(Hesse normal form)
    [acc,thetas,rhos]=hough_transform(A);
    %m is used to scale the representation of the hough space 
    m=max(acc(:));
    %resize the accumulator matrix from hough to the same size as image
    %only because it is easier for a human to look at. no other purpose
    [width,height]=size(A);
    resized=imresize(acc,[width,height]);
    %finds the peak in the accumulator matrix (for a human represented as the
    %bright spots on the accumulator image
    peaks=detect_max(acc);
    %creates a matrix same size as the image with the found lines drawn
    %into it
    line_matrix=create_line_matrix(peaks,rhos,thetas,A);
    %merges the image and the matrix lines together
    %the new line is not binary. so that the drawn lines can be drawn in a
    %differnt colour than the original edges. so they can be differentiated
    line_image=merge_lines_with_image(A,line_matrix);
    %results get shown, if the features are on the hough space image and the images with
    %the lines in it
    if(do_show_images)
        %Plot of the results
        f1 = figure;
        subplot(1,1,1), imshow(line_image)                       
            title('\fontsize{12}\bf Original image')
        hold on;
        %Plot of the hough space
        f2 = figure;
        subplot(1,1,1),imshow(resized,[0,m])                  
            title('\fontsize{12}\bf Hough image')
            xlabel('\fontsize{8} Winkel [°]')
            ylabel('\fontsize{8} Abstand [ ]')
    end
    %saves the image with the lines in the folder where the mat-file is.
    %same name as the original image except the line_ in front
    imwrite(line_image,strcat('lines_',image_name));
end

%creates the hough space (acc). out of the given binary image
%thetas is the angle vector in rad
%rhos is the distance vector in pixels
function [acc,thetas,rhos]=hough_transform(binary_image)
   
   %Are the edges white or black?
   isWhiteEdge= edge_value(binary_image);
   %rho and thta ranges
   thetas=degtorad(-90:90);
   [width,height]=size(binary_image);
   diag_len=round(sqrt(width*width+height*height));
   rhos=-diag_len:1:diag_len;
   
   %cache some reusable values
   cos_t = cos(thetas);
   sin_t=sin(thetas);
   num_thetas= size(thetas,2);
   
   %hough accumulator array of thetha vs rho
   acc=zeros(2*diag_len,num_thetas);
   %gets all the indexes from the pixels where an edge is
   if(isWhiteEdge)
    [y_idxs, x_idxs]= find(binary_image);
   else
    [y_idxs, x_idxs]= find(~binary_image);
   end
   
   %vote in hough accumulator
   for i= 1:size(x_idxs)
        x=x_idxs(i);
        y=y_idxs(i);
        for t_idx=1:num_thetas
            %calculate rho. diag_len is added for a positive index
            rho=round(x*cos_t(t_idx)+y*sin_t(t_idx))+diag_len;
            acc(rho,t_idx)=acc(rho,t_idx)+1;
        end
   end 
     
end

%gets the local max values and with that the lines
%return value is a number_peaks(max 10)x2. In this vectors there is rho and
%[rho,theta]
function [peaks]=detect_max(accumulator)

    %for detecting the local maximas in the accumulator matrix
    %an given matlab function was used. there is the name houghpeaks in it.
    %but could be used for any matrix to detect the local maximas. other
    %programming languages do have methods for that. not so matlab(or not
    %so for detecting them directly). We assumed that the main focus was to
    %lay on getting the hough space. so we used that method, especially
    %because of the hough transformation was one of the more challenging
    %algorithms to implement. 
    
    %in most already programmed hough algorithms there you can set a
    %threshold, value whatever to influence this local max detection. this
    %has e huge influence on which lines are recognised as lines and which
    %not. we decided the would go over the limits of this task
    peaks=houghpeaks(accumulator,10);
    
end

%this creates out of the hough space a matrix with the size of the image.
%in this matrix the found lines are drawn. the lines go trough the whole
%picture. it was not implemented to detect where the line would start an
%where it does end. 'cause in the "original" hough algorithm this is not
%given. nowadays this is normaly implemented -> but out of the scope of
%this task
function [line_matrix]=create_line_matrix(peaks,rhos,thetas,image)
    [width,height]=size(image);
    line_matrix=zeros(width,height);
    %draw those peaks on the image
    for element=1:size(peaks)
        peak=peaks(element,:);
        rho=rhos(peak(1));
        theta=thetas(peak(2));
        %get the beginn ing and endpoint of the line
        [x1,x2,y1,y2]=calculate_line_points(theta,rho,image);
        x = [x1 x2];   % x coordinates (running along matrix columns)
        y = [y1 y2];   % y coordinates (running along matrix rows)
        nPoints = round(sqrt(abs(diff(x))*abs(diff(x))+ abs(diff(y))*abs(diff(y))));  % Number of points in line
        rIndex = round(linspace(y(1), y(2), nPoints));          % Row indices
        cIndex = round(linspace(x(1), x(2), nPoints));          % Column indices       
        index = sub2ind(size(line_matrix), rIndex,cIndex);      % Linear indices
        line_matrix(index) = 255;  % Set the line pixels to the max value of 255 for uint8 types        
    end

end

%gives the start and enpoint of a straight line. start and endpoint of
%the edge of the image
%simply tries trough all 6 possible line directions. 
function [x1,x2,y1,y2]=calculate_line_points(theta,rho,image)
            
            [y,x]=size(image);
            theta=radtodeg(theta);
            if(theta==0)
                line([rho,rho],[1,y]);
            else
                %intersection x axis
                y1=1;
                x1=(sind(theta)/cosd(theta))*(rho/sind(theta) -y1);
                if(point_is_in_image(x1,y1,image))
                    %intersection y axis
                    x2=1;
                    y2=-cosd(theta)/sind(theta)*x2+rho/sind(theta);
                    if(~point_is_in_image(x2,y2,image))
                        %intersection y_parallel axis
                        x2=x;
                        y2=-cosd(theta)/sind(theta)*x2+rho/sind(theta);
                        if(~point_is_in_image(x2,y2,image))
                            %intersection x_parallel axis
                            y2=y;
                            x2=(sind(theta)/cosd(theta))*(rho/sind(theta) -y2);
                        end 
                    end
                else
                    %intersection y axis
                    x1=1;
                    y1=-cosd(theta)/sind(theta)*x1+rho/sind(theta);
                     if(point_is_in_image(x1,y1,image))
                        %intersection y_parallel axis
                        x2=x;
                        y2=-cosd(theta)/sind(theta)*x2+rho/sind(theta);
                        if(~point_is_in_image(x2,y2,image))
                            %intersection x_parallel axis
                            y2=y;
                            x2=(sind(theta)/cosd(theta))*(rho/sind(theta) -y2);
                        end 
                     else
                        %intersection y_parallel axis
                        x1=x;
                        y1=-cosd(theta)/sind(theta)*x1+rho/sind(theta);
                        if(point_is_in_image(x1,y1,image))
                            %intersection x_parallel axis
                            y2=y;
                            x2=(sind(theta)/cosd(theta))*(rho/sind(theta) -y2);    
                        end
                     end
                    
                end
            end
end

%returns 1 if the point is in the image 0 if it's not
function b=point_is_in_image(x,y,image)
        [y_max,x_max]=size(image);
       if(x>0&&x<=x_max)
           if(y>0&&y<=y_max)
               b=true;
           else
               b=false;
           end
       else
           b=false;
       end
end

%delets a possible frame around the image. if there is a frame it leads to
%false lines. returns the matrix without the frame
function [A] = delete_frame (A)
    [width,height]=size(A);
    %gets the value of the frame pixels
    frameValue=edge_value(A);
    %searches the upper left corner of the image without frame
    for row1 = 1:width
        for column1 = 1:height
            if (A(row1,column1) ~= frameValue)
                break
            end
        end
        if (A(row1,column1) ~= frameValue)
            break
        end
    end
    %searches the downer right corner of the image without frame
    for row2 = 0:width-1
        for column2 = 0:height-1
            if (A(width-row2,height-column2) ~= frameValue)
                break
            end
        end
        if (A(width-row2,height-column2) ~= frameValue)
            break
        end
    end
    %minimizes the image to the calculated inner corners
    A = A(row1:width-row2,column1:height-column2);        
end

%Calculates the searched value of images edges (balck-0/white-1):
function [value]= edge_value(A)
    numberOfzeros = sum(A(:) == 0);
    numberOfones = sum(A(:) == 1);
    if ( numberOfzeros >= numberOfones)
        value = 1;
    else
        value = 0;
    end
end

%returns a matrix in the image size with the detected lines in it drawn
function [line_image]=merge_lines_with_image(image,line_matrix)
    
     %white=255;
     %black=0;
     [n,m]=size(image);
     line_image=zeros(n,m,3);
     for row=1:n
         for column=1:m
             %if there is a detected line on this pixel. set the pixel to
             %green(best seen colour on white and black) so that it can be 
             %seen by eye
             if(line_matrix(row,column)==255)
                line_image(row,column,1)=0;
                line_image(row,column,2)=255;
                line_image(row,column,3)=0;
             else
                 %if there is no detected line set the original value of
                 %the image
                 line_image(row,column,1)=image(row,column);
                 line_image(row,column,2)=image(row,column);
                 line_image(row,column,3)=image(row,column);
             end
         end
     end
end


%% Erstelldatum: 09.04.2019
%  Verfasser: Delia De-Sassi, Patrick Wipfli