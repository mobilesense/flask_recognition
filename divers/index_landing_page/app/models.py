from app import app, db
import datetime

class Subscription(db.Model):
    __tablename__ = 'users'
    id = db.Column(db.Integer, primary_key=True) 
    date = db.Column(db.DateTime(timezone=False), default=datetime.datetime.utcnow, index=True)   
    mail = db.Column(db.String(120), unique=True)
            
    def __repr__(self):
        return '<Subscriptor %r>' % (self.username)

