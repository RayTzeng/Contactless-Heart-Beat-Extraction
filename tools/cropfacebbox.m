function [J,face,bbox,x,y] = cropfacebbox(img)
FaceDetect = vision.CascadeObjectDetector('FrontalFaceCART','MinSize',[150,150]);
bbox=FaceDetect.step(img);
if ~isempty(bbox)
    for i=1:size(bbox,1)
        J=imcrop(img,bbox(i,:));
    end
    face = 1;
    bboxPoints = bbox2points(bbox(1, :));
    x=mean(bboxPoints(:,1));
    y=mean(bboxPoints(:,2));
else
    J = img;
    face = 0;
    SZ=size(J);
    x=size(J,1)/2;
    y=size(J,2)/2;
end
end