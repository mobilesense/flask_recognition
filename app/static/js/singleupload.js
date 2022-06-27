// vars
var logo_holder = document.getElementById('single-holder'),
    icon_holder = document.getElementById('single-holder2'),
    logo_holder_img = document.getElementById('logo-holder-img'),
    icon_holder_img = document.getElementById('icon-holder-img'),    
    logo_holder_error = "holder-error",
    icon_holder_error = "holder-error2",
    tests = {
      filereader: typeof FileReader != 'undefined',
      dnd: 'draggable' in document.createElement('span'),
      formdata: !!window.FormData
    }, 
    /*support = {
      filereader: document.getElementById('single-filereader')
    },*/
    acceptedTypes = {'image/png': true, 'image/jpeg': true, 'image/jpg': true, 'image/gif': true},
    logo_fileupload = document.getElementById('single-upload'),
    icon_fileupload = document.getElementById('single-upload2');
// end vars

function previewfile(file, holder, holder_img, holder_error) {
    if (tests.filereader === true && acceptedTypes[file.type] === true) {
        var reader = new FileReader();
        reader.onload = function (event) {
            var image = new Image();
            image.src = event.target.result;
            //image.width = 250; // a fake resize
            holder_img.innerHTML = "";
            holder_img.appendChild(image);
            return image;
        };
        reader.readAsDataURL(file);
        
        //send file or reader ?
        console.log("sending..");        
    }
    else {
        // error on file type
        console.log("error format");        
        $(holder_error).text('Not accepted file type. Have to be: png, jpeg, jpg, or gif.');
    }

}

function upload(holder, holder_img, holder_error, fileupload){
    if (tests.dnd) { 
      holder.ondragover = function () { this.className = 'hover dropzone '; return false; };
      holder.ondragend = function () { 
        this.className = ''; 
        return false;         
      };      
      holder.ondrop = function (e) {
        //this.className = '';
        e.preventDefault();
        previewfile(e.dataTransfer.files[0], holder, holder_img, holder_error); 
      }
    } else {
      fileupload.className = 'hidden';
      fileupload.querySelector('input').onchange = function () {
        previewfile(this.files[0], holder, holder_img, holder_error);
      };
    }
}

upload(logo_holder, logo_holder_img, logo_holder_error, logo_fileupload);
upload(icon_holder, icon_holder_img, icon_holder_error, icon_fileupload);



