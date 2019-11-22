## Installation of FunGAPsite-webapp
Please install FunGAP first, if you haven't install, click [here](https://github.com/CompSynBioLab-KoreaUniv/FunGAP).

Download FunGAPsite-webapp from github
```
cd $FUNGAP_DIR #the path where FunGAP is installed
git clone https://github.com/takimailto/FunGAPsite-webapp.git
cd FunGAPsite-webapp
conda activate fungap #The virtual environment you used to install FunGAP
pip install -r requirements.txt
```
#### Download and install database mysql
```
sudo apt-get update
sudo apt-get install mysql-server mysql-client
```
Set the root user and password. Then, run mysql
```
mysql -u root -p
```
Type password, create a database to store message from the FunGAPsite-webapp.
For example:
```
CREATE DATABASE FunGAPsite CHARACTER SET utf8
```
You can type ```SHOW DATABASES;``` to see the existing database.
Then, change the django setting.
```
cd FunGAPsite
vi settings.py
```
change
```DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
        'NAME': '', # name of the database you create
        'HOST': '', # 127.0.0.1 by default
        'USER': '', # root
        'PASSWORD': '', # the password you set
        'PORT': '', # 3306 by defult
    }
}
```
Save the change. 
#### Change the configuration
Go to the ```fungap``` directory
```
cd ..
cd fungap
```
Change ```views.py```:
```
BASE_DIR = '***/FunGAPsite-webapp/FunGAPsite/work_space/users'  #***: The FunGAP path in your computer 
```
Change ```tasks.py```:
```
FunGAP_dir = '***'  # ***: The FunGAP path in your computer
```
Make migration to your database by run:
```
python manage.py makemigrations
python manage.py migrate
```
Now, you can run the server by:
```
python manage.py runserver
```
Open the URL as prompted， you can start your project now.

#### Create superuser
You can also create superuser for more convenient management by:
```
python manage.py createsuperuser
```
Then follow the instructions.
