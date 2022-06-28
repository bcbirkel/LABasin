clc, clf, clear
%dungeon.m
%This is a mini role-playing game that cannot be saved.
%You first need to enter a game number that is a positive integer.
%Then you need to input your name and race.  
%Your race affects your stats so choose wisely.  Your default stats range
%from 3 to 18.
%
%   Humans have default stats.  They get no bonuses.
%
%   Elves have lower strength and stamina, but greater dexterity and
%   intelligence.
%
%   Dwarves have lower dexterity and intelligence, but greater strength
%   and stamina.
%
%   Halflings have greater dexterity and wisdom, but less strength and
%   stamina.
%
%   There is a secret race, so you can play around with various numbers and
%   names to try and unlock it.
%
%When inputting infomation please be sure to use only lower case charcters.
%
%Once you are happy with your stats, then you must choose a class and deity.  
%The class and deity affect what skills you have and how you level up.  
%
%   Warriors have greater attacking power and greater defensive power.
%   They can equip the best regular weapons and armor.  They have stronger
%   healing magic and only get one spell.
%
%   Rogues are the only class that can initially steal.  They have decent hit
%   percentages too.  They get only one spell and have stronger attacking
%   magic power.
%
%   Wizards have the strongest attack magic.  Their physical prowess
%   however is somewhat lacking.  They get two spells, at least one must
%   be attacking.
%
%   Priests have the strongest healing magic.  They are slightly stronger
%   than wizards physically, but they are the most challenging class to
%   play.  They get two spells, at least one must be healing.
%
%   There is also a hidden class.  It is the strongest class.  You can get
%   it if you input the proper name or number at the prompt.
%
%The five gods are of water, fire, thunder, ice, and earth.  
%Then you must choose your original spell or spells.
%
%Next your character goes shopping for basic equipment.  Different weapons
%have advantages against different monsters so the most expensive one is not
%always the best, but generally the more expensive one does have a stronger attack.
%Armor always increases strength with cost.  Be careful because you can
%only have one set of armor and one weapon at a time.
%
%Then the adventure begins.
%Use the numberpad to move around the map when prompted.
%At the same prompt you can enter a,s,z,q, or w  to check a variety of
%character information.
%
%   a checks current inventory including gold, healing herbs, magic
%   tonics, weapon, and armor.  It also shows the strength of your weapon
%   and armor, which increases as you increase your level.
%
%   s shows your status.  This includes current health and magic power as
%   well as level and experience.  Every 100 points of experience
%   increases your level increasing your stats.
%
%   z gives your fundamental stats along with race, class, and deity.
%
%   q shows your magic capabilities and each spell's cost.  As you increase
%   your level you will be given the option to increase your magic options
%   or gain the rogue's ability to steal.
%
%   w shows your battle bonuses which are based on class and basic stats.
%   They increase with your level.
%
%Use the e and f keys to perform healing tasks while outside of battle.   
%
%   e will use an item. 
%   
%   f will use cure magic if available.
%
%   Use 0 to quit.
%
%Battles and other events happen randomly.  These events will gain you
%experience to increase your level so you are strong enough to destroy the
%evil Devousius.  Should you die or defeat Devousius it wil end the
%program.

help dungeon
go_on=input('Hit enter when ready to continue','s');
clc

%This game number section resets the random list.
game_number=input('Please input a game number.  It may be any positive integer.   ')
if game_number>0
    game_number=ceil(game_number);
    rand(game_number);
else
    game_number=1;
    rand(game_number);
end
clc

%This section welcomes the player and creates some basic information
%name
disp('Welcome to the dungeon of Devousius')
name=input('What is your name, mighty hero?    ','s');

%race
disp('The races of this world are humans, elves, dwarves, and halflings.')
race_flag=0;
humans=0;
elves=1;
dwarves=2;
halflings=3;
half_elf=45;
fprintf('Welcome %s, live with honor.\n',name)

while race_flag==0
race=input('Who are your closest kin?  0- humans; 1- elves; 2-dwarves; 3-halflings   ');
c=ceil(20*rand(1));
for d=1:c
   r=rand(d);
end
%This block generates the five main stats to be used throughout the game
strength=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));
stamina=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));
dexterity=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));
wisdom=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));
intelligence=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));

%Each race has different stat changes defined below
if race==1
    dexterity=dexterity+1;
    intelligence=intelligence+1;
    stamina=stamina-1;
    strength=strength-1;
    race_flag=1;
elseif race==0
    race_flag=1;
elseif race==2
    strength=strength+1;
    stamina=stamina+1;
    dexterity=dexterity-1;
    intelligence=intelligence-1;
    race_flag=1;
elseif race==3
    strength=strength-1;
    dexterity=dexterity+1;
    wisdom=wisdom+1;
    stamina=stamina-1;
    race_flag=1;
elseif race==45
    dexterity=dexterity+1;
    intelligence=intelligence+1;
    race_flag=1;
else
    race=floor(4*rand(1));
    race_flag=1;
end
end

%This displays the stats and asks if it should get new stats
fprintf('Your strength is %1.0f.\nYour stamina is %1.0f.\nYour dexterity is %1.0f.\nYour intelligence is %1.0f.\nYour wisdom is %1.0f.\n',strength,stamina,dexterity,intelligence,wisdom')
no=0;
n=0;
yes=1;
y=1;
reroll=input('Are these acceptable values?  ');
clc
while reroll==0
    strength=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));
    stamina=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));
    dexterity=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));
    wisdom=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));
    intelligence=ceil(6*rand(1))+ceil(6*rand(1))+ceil(6*rand(1));
    clc
    %Each race has different stat changes defined below
if race==1
    dexterity=dexterity+1;
    intelligence=intelligence+1;
    stamina=stamina-1;
    strength=strength-1;
elseif race==2
    strength=strength+1;
    stamina=stamina+1;
    dexterity=dexterity-1;
    intelligence=intelligence-1;
elseif race==3
    strength=strength-1;
    dexterity=dexterity+1;
    wisdom=wisdom+1;
    stamina=stamina-1;
elseif race==45
    dexterity=dexterity+1;
    intelligence=intelligence+1;
end
    fprintf('Your strength is %1.0f.\nYour stamina is %1.0f.\nYour dexterity is %1.0f.\nYour intelligence is %1.0f.\nYour wisdom is %1.0f.\n\n',strength,stamina,dexterity,intelligence,wisdom')
    reroll=input('Are these acceptable values?  ');
end

%class
disp('The classes available are warrior, rogue, wizard, or priest.')
warrior=0;
rogue=1;
wizard=2;
priest=3;
holy_knight=27;
class_flag=0;
while class_flag==0
class=input('Under which discipline do you train?  0- warrior; 1- rogue; 2- wizard; 3- priest   ');
if class==0 | class==1 | class==2 | class==3 | class==27
    class_flag=1;
else
   if strength>=dexterity & strength>=intelligence & strength>=wisdom
       class=0;
       class_flag=1;
   elseif dexterity>=strength & dexterity>=intelligence & dexterity>=wisdom
       class=1;
       class_flag=1;
   elseif intelligence>=strength & intelligence>=dexterity & intelligence>=wisdom
       class=2;
       class_flag=1;
   elseif wisdom>=strength & wisdom>=dexterity & wisdom>=intelligence
       class=3;
       class_flag=1;
   end
end
end
%deity
disp('The deities of this world are the serpent Leviathon, the bird Pheonix,')
disp('the elder kin Ramuh, the giant Titan, and the warrior Shiva,')
leviathon=0;
pheonix=1;
ramuh=2;
titan=3;
shiva=4;
deity_flag=0;
c=ceil(20*rand(1));
    for d=1:c
        r=rand(d);
    end
while deity_flag==0
deity=input('Under whose star do you pray? 0- Leviathon; 1- Pheonix; 2- Ramuh; 3- Titan; 4-Shiva   ');
clc
if deity==0 | deity==1 | deity==2 | deity==3 | deity==4
    deity_flag=1;
else
    deity=floor(5*rand(1));
    deity_flag=1;
end
end
%Now we need to define the character's starting stats.  These are based on
%class and stats.
%First for warrior
if class==0
    health=stamina+ceil(6*rand(1));
    magic=ceil(wisdom/3+ceil(2*rand(1)));
    attackbonus=ceil(strength/6);
    hitbonus=ceil(dexterity/6);
    magebonus=ceil(intelligence/12);
    curebonus=ceil(wisdom/6);
    armorclass=2;
    weaponclass=3;
    steal_option=0;
%Next for rogue
elseif class==1
    health=ceil(stamina/2+ceil(8*rand(1)));
    magic=ceil(wisdom/2+ceil(2*rand(1)));
    attackbonus=ceil(strength/8);
    hitbonus=ceil(dexterity/4);
    magebonus=ceil(intelligence/6);
    curebonus=ceil(wisdom/12);
    armorclass=2;
    weaponclass=2;
    steal_option=1;
%Next up is wizard
elseif class==2
    health=ceil(stamina/4+ceil(4*rand(1)));
    magic=ceil(wisdom+ceil(4*rand(1)));
    attackbonus=ceil(strength/12);
    hitbonus=ceil(dexterity/12);
    magebonus=ceil(intelligence/3);
    curebonus=ceil(wisdom/4);
    armorclass=0;
    weaponclass=1;
    steal_option=0;
%Next is priest
elseif class==3
    health=ceil(stamina/3+ceil(5*rand(1)));
    magic=ceil(wisdom+ceil(4*rand(1)));
    attackbonus=ceil(strength/8);
    hitbonus=ceil(dexterity/8);
    magebonus=ceil(intelligence/4);
    curebonus=ceil(wisdom/3);
    armorclass=0;
    weaponclass=1;
    steal_option=0;
%Last is the holy knight
elseif class==27
    health=ceil(stamina+ceil(8*rand(1)));
    magic=ceil(wisdom+ceil(6*rand(1)));
    attackbonus=ceil(strength/5);
    hitbonus=ceil(dexterity/5);
    magebonus=ceil(intelligence/3);
    curebonus=ceil(wisdom/3);
    armorclass=3;
    weaponclass=3;
    steal_option=0;
end

%Now we need to put in magic
%These defined variables will be used later and changed in the blocks
healing=0;
attacking=1;
cure_option=0;
attack_option=0;
cure_option_strong=0;
attack_option_strong=0;
level_skill=0;

%Warrior gets one spell
if class==0
    disp('Your god will give you one spell to keep you safe.')
    spell_1=input('Which do you want? 0- healing; 1- attacking  ');
    if spell_1==0
        disp('You now have the basic cure spell. Use it wisely.')
        cure_option=1;
    elseif spell_1==1
        attack_option=1;
        if deity==0
            disp('You now have the water spell.  Use it wisely.')
        elseif deity==1
            disp('You now have the fire spell.  Use it wisely.')
        elseif deity==2
            disp('You now have the thunder spell.  Use it wisely.')
        elseif deity==3
            disp('You now have the earth spell.  Use it wisely.')
        elseif deity==4
            disp('You now have the ice spell.  Use it wisely.')
        end
    else
        disp('You now have the basic cure spell. Use it wisely.')
        cure_option=1;
    end
end

%Rogue gets one spell
if class==1
    disp('Your god will give you one spell to keep you safe.')
    spell_1=input('Which do you want? 0- healing; 1- attacking  ');
    if spell_1==0
        disp('You now have the basic cure spell. Use it wisely.')
        cure_option=1;
    elseif spell_1==1
        attack_option=1;
        if deity==0
            disp('You now have the water spell.  Use it wisely.')
        elseif deity==1
            disp('You now have the fire spell.  Use it wisely.')
        elseif deity==2
            disp('You now have the thunder spell.  Use it wisely.')
        elseif deity==3
            disp('You now have the earth spell.  Use it wisely.')
        elseif deity==4
            disp('You now have the ice spell.  Use it wisely.')
        end
    else
        attack_option=1;
        if deity==0
            disp('You now have the water spell.  Use it wisely.')
        elseif deity==1
            disp('You now have the fire spell.  Use it wisely.')
        elseif deity==2
            disp('You now have the thunder spell.  Use it wisely.')
        elseif deity==3
            disp('You now have the earth spell.  Use it wisely.')
        elseif deity==4
            disp('You now have the ice spell.  Use it wisely.')
        end
    end
end
%Wizards get two spells, at least one is offensive
if class==2
    disp('Your god has bestowed upon you two spells to keep you safe.')
    spell_1=input('Which do you want? 0- healing; 1- attacking  ');
    if spell_1==0
        disp('You now have the basic cure spell. Use it wisely.')
        disp('Your second spell must be attacking.')
        cure_option=1;
        attack_option=1;
        if deity==0
            disp('You now have the basic water spell.  Use it wisely.')
        elseif deity==1
            disp('You now have the basic fire spell.  Use it wisely.')
        elseif deity==2
            disp('You now have the basic thunder spell.  Use it wisely.')
        elseif deity==3
            disp('You now have the basic earth spell.  Use it wisely.')
        elseif deity==4
            disp('You now have the basic ice spell.  Use it wisely.')
        end
    elseif spell_1==1
        attack_option=1;
        if deity==0
            disp('You now have the basic water spell.  Use it wisely.')
        elseif deity==1
            disp('You now have the basic fire spell.  Use it wisely.')
        elseif deity==2
            disp('You now have the basic thunder spell.  Use it wisely.')
        elseif deity==3
            disp('You now have the basic earth spell.  Use it wisely.')
        elseif deity==4
            disp('You now have the basic ice spell.  Use it wisely.')
        end
        spell_2=input('What will be your second spell? 0- healing; 1- attacking  ');
        if spell_2==1
            attack_option_strong=1;
            if deity==0
            disp('You now have a stronger water spell.  Use it wisely.')
            elseif deity==1
            disp('You now have a stronger fire spell.  Use it wisely.')
            elseif deity==2
            disp('You now have a stronger thunder spell.  Use it wisely.')
            elseif deity==3
            disp('You now have a stronger earth spell.  Use it wisely.')
            elseif deity==4
            disp('You now have a stronger ice spell.  Use it wisely.')
            end   
        elseif spell_2==0
               disp('You now have the basic cure spell. Use it wisely.')
               cure_option=1;
        else
               attack_option_strong=1;
            if deity==0
            disp('You now have a stronger water spell.  Use it wisely.')
            elseif deity==1
            disp('You now have a stronger fire spell.  Use it wisely.')
            elseif deity==2
            disp('You now have a stronger thunder spell.  Use it wisely.')
            elseif deity==3
            disp('You now have a stronger earth spell.  Use it wisely.')
            elseif deity==4
            disp('You now have a stronger ice spell.  Use it wisely.')
            end   
        end
        else
    attack_option=1;
        if deity==0
            disp('You now have the basic water spell.  Use it wisely.')
        elseif deity==1
            disp('You now have the basic fire spell.  Use it wisely.')
        elseif deity==2
            disp('You now have the basic thunder spell.  Use it wisely.')
        elseif deity==3
            disp('You now have the basic earth spell.  Use it wisely.')
        elseif deity==4
            disp('You now have the basic ice spell.  Use it wisely.')
        end
           attack_option_strong=1;
            if deity==0
            disp('You now have a stronger water spell.  Use it wisely.')
            elseif deity==1
            disp('You now have a stronger fire spell.  Use it wisely.')
            elseif deity==2
            disp('You now have a stronger thunder spell.  Use it wisely.')
            elseif deity==3
            disp('You now have a stronger earth spell.  Use it wisely.')
            elseif deity==4
            disp('You now have a stronger ice spell.  Use it wisely.')
            end   
     end
end

%priests get two spells at least one healing
if class==3
    disp('Your god has bestowed upon you two spells to keep you safe.')
    spell_1=input('Which do you want? 0- healing; 1- attacking  ');
    if spell_1==0
        disp('You now have the basic cure spell. Use it wisely.')
        cure_option=1;
        spell_2=input('What will be your second spell? 0- healing; 1- attacking  ');
        if spell_2==0
            disp('You now have a stronger cure spell.  Use it wisely.')
            cure_option_strong=1;
        elseif spell_2==1
            attack_option=1;
            if deity==0
            disp('You now have the basic water spell.  Use it wisely.')
            elseif deity==1
            disp('You now have the basic fire spell.  Use it wisely.')
            elseif deity==2
            disp('You now have the basic thunder spell.  Use it wisely.')
            elseif deity==3
            disp('You now have the basic earth spell.  Use it wisely.')
            elseif deity==4
            disp('You now have the basic ice spell.  Use it wisely.')
            end
        else
            disp('You now have a stronger cure spell.  Use it wisely.')
            cure_option_strong=1;
        end
    elseif spell_1==1
        attack_option=1;
        if deity==0
            disp('You now have the basic water spell.  Use it wisely.')
        elseif deity==1
            disp('You now have the basic fire spell.  Use it wisely.')
        elseif deity==2
            disp('You now have the basic thunder spell.  Use it wisely.')
        elseif deity==3
            disp('You now have the basic earth spell.  Use it wisely.')
        elseif deity==4
            disp('You now have the basic ice spell.  Use it wisely.')
        end
        disp('Your second spell must be healing.')
        disp('You now have the basic cure spell. Use it wisely.')
        cure_option=1;
    else
        disp('You now have the basic cure spell. Use it wisely.')
        cure_option=1;
        disp('You now have a stronger cure spell.  Use it wisely.')
        cure_option_strong=1;
    end
end
%Holy knight gets two healing spell and two attacking spells
if class==27
    disp('You now have the basic cure spell. Use it wisely.')
    cure_option=1;
    disp('You now have a stronger cure spell.  Use it wisely.')
    cure_option_strong=1;
    attack_option=1;
        if deity==0
            disp('You now have the basic water spell.  Use it wisely.')
        elseif deity==1
            disp('You now have the basic fire spell.  Use it wisely.')
        elseif deity==2
            disp('You now have the basic thunder spell.  Use it wisely.')
        elseif deity==3
            disp('You now have the basic earth spell.  Use it wisely.')
        elseif deity==4
            disp('You now have the basic ice spell.  Use it wisely.')
        end
    attack_option_strong=1;
        if deity==0
            disp('You now have a stronger water spell.  Use it wisely.')
        elseif deity==1
            disp('You now have a stronger fire spell.  Use it wisely.')
        elseif deity==2
            disp('You now have a stronger thunder spell.  Use it wisely.')
        elseif deity==3
            disp('You now have a stronger earth spell.  Use it wisely.')
        elseif deity==4
            disp('You now have a stronger ice spell.  Use it wisely.')
        end   
    end
%Now a pause so to display the text.
go_on=input('Hit enter when you are ready to continue.','s'  );
clc

%We set starting values for experience and gold.
level=1;
experience=0;
gold=ceil(20*rand(1)+20*rand(1)+20*rand(1)+20*rand(1)+20*rand(1)+20*rand(1));

%Now we show the starting stats
disp('Here is your character')
disp(name)
if class==0
    disp('Warrior')
elseif class==1
    disp('Rogue')
elseif class==2
    disp('Wizard')
elseif class==3
    disp('Priest')
elseif class==27
    disp('Holy Knight')
end
if race==0
    disp('Human')
elseif race==1
    disp('Elf')
elseif race==2
    disp('Dwarf')
elseif race==3
    disp('Halfling')
elseif race==45
    disp('Half-elf')
end
go_on4=input('Hit enter when ready to continue','s')
clc
%show stats
fprintf('\nYour strength is %1.0f.\nYour stamina is %1.0f.\nYour dexterity is %1.0f.\nYour intelligence is %1.0f.\nYour wisdom is %1.0f.\n\n',strength,stamina,dexterity,intelligence,wisdom')
go_on3=input('Hit enter when ready to continue','s')
clc
%Show other stuff
fprintf('\nYour starting max health is %1.0f.\nYour starting max magic points is %1.0f.\nYour attack damage bonus is %1.0f.\nYour to hit bonus is %1.0f.\nYour magic attack bonus is %1.0f.\nYour cure magic bonus is %1.0f.\n',health,magic,attackbonus,hitbonus,magebonus,curebonus)
fprintf('\nYou are starting at level %1.0f with %1.0f experience.\nYou have %1.0f gold to spend.\n',level,experience,gold)
go_on2=input('Hit enter when ready to continue.','s')
clc

%Now the character needs to buy some basics.
fprintf('You have %1.0f gold, time to go shopping.\n', gold)
armory=0;
weaponry=1;
item_shop=2;
shopping=0;
healing_herbs=0;
magic_tonic=0;
armor_type=3;
weapon_name=0;
weapontype=1;

while shopping~=3
    shopping=input('Where do you want to go?  0- Armory; 1- Weaponry; 2- Item Shop?  Hit enter when done. ');
    %buy some armor
    if shopping==0
        disp('We have various armors available, which would you like?')
        %warrior and holy knight
        if class==0 | class==27
        armor=input('1- leather armor 10 gold; 2- chain mail 30 gold; 3- plate mail 80 gold   ');
            if armor==1
                if gold >=10
                gold=gold-10
                armorclass=3;
                armor_type=0;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==2
                if gold >=30
                gold=gold-30
                armorclass=5;
                armor_type=1;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==3
                if gold >=80
                gold=gold-80
                armorclass=9;
                armor_type=2;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
                gold=gold
                armmorclass=armorclass;
                armor_type=armor_type;
            end
            %rogue
        elseif class==1
        armor=input('1- cape 5 gold; 2- leather armor 10 gold; 3- light chain mail 40 gold   ');
            if armor==1
                if gold >=5
            gold=gold-5
            armorclass=3;
            armor_type=4;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==2
                if gold >=10
            gold=gold-10
            armorclass=5;
            armor_type=0;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==3
                if gold >=40
            gold=gold-40
            armorclass=7;
            armor_type=5;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
            gold=gold
            armorclass=armorclass;
            armor_type=armor_type;
            end
            %wizard and priest
        elseif class==2 | class==3
            armor=input('1-robe 5 gold; 2- cloak 15 gold; 3- padded robe 25 gold   ');
             if armor==1
                 if gold >=5
            gold=gold-5
            armorclass=1;
            armor_type=6;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==2
                if gold >=15
            gold=gold-15
            armorclass=2;
            armor_type=7;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==3
                if gold >=25
            gold=gold-25
            armorclass=4;
            armor_type=8;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
            gold=gold
            armorclass=armorclass;
            armor_type=armor_type;
            end
        end %armors
        %weapons
    elseif shopping==1
        disp('We have weapons for every discipline, take your pick.')
        if class==0 | class==27
            weapon=input('1- long sword 25 gold; 2- mace 30 gold; 3- spear 40 gold   ');
            if weapon==1
                if gold >=25
                gold=gold-25
                weaponclass=6;
                %weapontype is slashing (0), bashing(1) or piercing(2)
                weapontype=0;
                weapon_name=1;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==2
                if gold >=30
                gold=gold-30
                weaponclass=6;
                weapontype=1;
                weapon_name=2;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==3
                if gold >=40
                gold=gold-40
                weaponclass=7;
                weapontype=2;
                weapon_name=3;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
                gold=gold
                weaponclass=3;
            end
        elseif class==1
            weapon=input('1- knife 10 gold; 2- club 15 gold; 3- rapier 30 gold   ');
            if weapon==1
                if gold >=10
                gold=gold-10
                weaponclass=4;
                weapontype=0;
                weapon_name=4;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==2
                if gold >=15
                gold=gold-15
                weaponclass=4;
                weapontype=1;
                weapon_name=5;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==3   
                if gold >=30
                gold=gold-30
                weaponclass=6;
                weapontype=2;
                weapon_name=6;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
                gold=gold
                weaponclass=2;
            end
        elseif class==2 | class==3
             weapon=input('1- knife 10 gold; 2- staff 20 gold; 3- dagger 15 gold   ');
            if weapon==1
                if gold >=10
                gold=gold-10
                weaponclass=4;
                weapontype=0;
                weapon_name=4;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==2
                if gold >=20
                    gold=gold-20
                    weaponclass=5;
                    weapontype=1;
                    weapon_name=7;
                    else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==3
                if gold >=15
                    gold=gold-15
                    weaponclass=4;
                    weapontype=2;
                    weapon_name=8;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
                gold=gold
                weaponclass=1;
            end
        end %weapons
        %items
        elseif shopping==2
            disp('Welcome, we have healing herbs and magic tonics, how can we help you?')
            item=input('1- healing herb 10 gold; 2- magic tonic 50 gold  ');
            if item==1
                if gold >= 10
                    gold=gold-10
                    disp('Thank you for buying a healing herb.')
                    healing_herbs=healing_herbs+1;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif item==2
                if gold >=50
                    disp('Thank you for buying a magic tonic.')
                    gold=gold-50
                    magic_tonic=magic_tonic+1;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            end
        end
    end
    clc
    %shopping's done
    %Now we use a figues as a guide for the character
%First we define some starting points
x_home=5;
y_home=5;
x_dungeon=12;
y_dungeon=17;
x_mall=15;
y_mall=3;
x_position=10;
y_position=10;
exit_location_x=ceil(20*rand(1));
exit_location_y=ceil(20*rand(1));
boss_location_x=ceil(20*rand(1))+ceil(20*rand(1))-1;
boss_location_y=ceil(20*rand(1))+ceil(20*rand(1))-1;
treasure_location_1x=ceil(39*rand(1));
treasure_location_1y=ceil(39*rand(1));
treasure_location_2x=ceil(39*rand(1));
treasure_location_2y=ceil(39*rand(1));
x_position2=exit_location_x;
y_position2=exit_location_y+1;
playing_game=1;
holy_armor=0;
c=ceil(20*rand(1));
    for d=1:c
        r=rand(d);
    end
    
%Then we define some variables for later
motion=5;
current_health=health;
current_magic=magic;
location=1;
s=10;
a=11;
z=12;
w=13;
q=14;
e=15;
f=16;
orcs_action=5;
boss_health=150;
plot(x_home,y_home,'s',x_dungeon,y_dungeon,'p',x_mall,y_mall,'*',x_position,y_position,'d')
legend('Home','Dungeon','Mall',name)
axis([0 20 0 20])
grid on
    
while current_health>0 & boss_health>0 & playing_game==1
c=ceil(20*rand(1));
    for d=1:c
        r=rand(d);
    end
    %Now we change the location of the character
if location==1 %begin while
    motion=input('Use the numberpad to move around the map.   ');
    random_event=round(20*rand(1));
    if motion==1  %begin motion
        x_position=x_position-1;
        y_position=y_position-1;
    elseif motion==2
          y_position=y_position-1;
          x_position=x_position;
    elseif motion==3
          y_position=y_position-1;
          x_position=x_position+1;
     elseif motion==4
          y_position=y_position;
          x_position=x_position-1;
     elseif motion==5
          y_position=y_position;
          x_position=x_position;
     elseif motion==6
          y_position=y_position;
          x_position=x_position+1;
     elseif motion==7
          y_position=y_position+1;
          x_position=x_position-1;
     elseif motion==8
          y_position=y_position+1;
          x_position=x_position;
     elseif motion==9
          y_position=y_position+1;
          x_position=x_position+1;
      elseif motion==0
          y_position=y_position;
          x_position=x_position;
          clf,clc
          random_event=-12;
          playing_game=0;
          disp('Thanks for playing, good-bye')
      elseif motion==16
          if cure_option==0
             disp('You do not know that spell.')
                               
          elseif cure_option==1 & cure_option_strong==0 & current_magic>=3
                                disp('You feel much better.')
                                current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                current_magic=current_magic-3;
                                
                                if current_health>health
                                    current_health=health;
                                end
                            
                                fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                
          elseif cure_option==1 & cure_option_strong==0 & current_magic<3
                                disp('You do not have enough magic!')
                                
          elseif cure_option_strong==1
                                cure_strength=input('Which cure will you use? 1- weak cure or 2- strong cure  ');
                                
                                if cure_strength==1 & current_magic>=3
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                    current_magic=current_magic-3;
                                
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                 elseif cure_strength==1 & current_magic<3
                                    disp('You do not have enough magic!')
                                    
                                 elseif cure_strength==2 & current_magic>=5
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(6*rand(1))+ceil(level*rand(1));
                                    current_magic=current_magic-5;
                                 
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==2 & current_magic<5
                                   disp('You do not have enough magic!')
                                 end  %strong cure magic
          end%cure magic
      elseif motion==15
          fprintf('You have %1.0f healing herbs and %1.0f magic tonics.\n',healing_herbs,magic_tonic)
                            item_action=input('What will you use? 1- healing herb; 2- magic tonic?   ');

                            if item_action==1 & healing_herbs>=1                             
                                healing_herbs=healing_herbs-1;
                                disp('You feel better')
                                current_health=current_health+ceil(10*rand(1))+curebonus;
                                
                                if current_health>health
                                  current_health=health;
                                end
                                
                            elseif item_action==1 & healing_herbs<1
                                disp('You do not have any left.')
                                
                            elseif item_action==2 & magic_tonic>=1
                                magic_tonic=magic_tonic-1;
                                current_magic=current_magic+ceil(10*rand(1))+curebonus+magebonus;
                                disp('You feel energized')
                                if current_magic>magic
                                    current_magic=magic;
                                end
                            elseif item_action==2 & magic_tonic<1
                                disp('You do not have any left.')
                            end  %items
      elseif motion==10
          fprintf('\n%s\nYour current health is %1.0f/%1.0f.\nYour current magic power is %1.0f/%1.0f.\nYour current gold is %1.0f.\nYour current level is %1.0f\nYour current experience is %1.0f.\n',name,current_health,health,current_magic,magic,gold,level,experience)
      elseif motion==11
          fprintf('You currently have %1.0f gold.\nYou have %1.0f healing herbs.\nYou have %1.0f magic tonics.\n',gold,healing_herbs,magic_tonic)
          if weapon_name==0
              disp('You are fighting bare-handed.')
          elseif weapon_name==1
              disp('You have a long sword.')
          elseif weapon_name==2
              disp('You weapon is a mace.')
          elseif weapon_name==3
              disp('A long spear is your weapon.')
          elseif weapon_name==4
              disp('Your weapon is a knife.')
          elseif weapon_name==5
              disp('You have a club.')
          elseif weapon_name==6
              disp('Your weapon of choice is a rapier.')
          elseif weapon_name==7
              disp('You are carrying a staff.')
          elseif weapon_name==8
              disp('You are holding a dagger.')
          elseif weapon_name==9
              disp('You are weilding the Demonslayer.')
          end
          if armor_type==3
              disp('You have no armor.')
          elseif armor_type==0
              disp('You have leather armor.')
          elseif armor_type==1
              disp('You have chain mail.')
          elseif armor_type==2
              disp('You are donning plate mail.')
          elseif armor_type==4
              disp('You are wearing a cape.')
          elseif armor_type==5
              disp('You have light chain mail.')
          elseif armor_type==6
              disp('You are wearing a robe.')
          elseif armor_type==7
              disp('Your armor is a cloak.')
          elseif armor_type==8
              disp('You have a padded robe.')
          elseif armor_type==9
              disp('The Holy Armor is protecting you.')
          end
          fprintf('\nThe strength of you weapon is %1.0f.\nThe strength of your armor is %1.0f.\n\n',weaponclass,armorclass)
      elseif motion==12
          fprintf('Your strength is %1.0f.\nYour stamina is %1.0f.\nYour dexterity is %1.0f.\nYour intelligence is %1.0f.\nYour wisdom is %1.0f.\n',strength,stamina,dexterity,intelligence,wisdom)
          if race==0
              disp('You are a human.')
          elseif race==1
              disp('You are an elf.')
          elseif race==2
              disp('You are a dwarf.')
          elseif race==3
              disp('You are a halfling.')
          elseif race==45
              disp('You are a half-elf.')
          end
          if class==0
              disp('You are a warrior.')
          elseif class==1
              disp('You are a rogue.')
          elseif class==2
              disp('You are a wizard.')
          elseif class==3
              disp('You are a priest.')
          elseif class==27
              disp('You are a holy knight.')
          end
          
          if deity==0
              disp('Your guardian is the serpent Leviathon.')
          elseif deity==1
              disp('Your life is the ash of the Pheonix.')
          elseif deity==2
              disp('The power of Ramuh flows through you.')
          elseif deity==3
              disp('Your strength is that of Titan.')
          elseif deity==4
              disp('Shiva blesses you.')
          end
      elseif motion==13
          fprintf('Your attack damage bonus is %1.0f.\nYour to hit bonus is %1.0f.\nYour magic attack bonus is %1.0f.\nYour cure magic bonus is %1.0f.\n',attackbonus,hitbonus,magebonus,curebonus)
      elseif motion==14
          if cure_option==1
              disp('You have a basic cure spell. 3 magic per use.')
          end
          if cure_option_strong==1
              disp('You have a strong cure spell. 5 magic per use.')
          end
          if attack_option==1
              if deity==0
                  disp('You have a basic water spell. 4 magic per use.')
              elseif deity==1
                  disp('You have a basic fire spell. 4 magic per use.')
              elseif deity==2
                  disp('You have a basic lightning spell. 4 magic per use.')
              elseif deity==3
                  disp('You have a basic earth spell. 4 magic per use.')
              elseif deity==4
                  disp('You have a basic ice spell. 4 magic per use.')
              end
          end
          if attack_option_strong==1
              if deity==0
                  disp('You have a strong water spell. 7 magic per use.')
              elseif deity==1
                  disp('You have a strong fire spell. 7 magic per use.')
              elseif deity==2
                  disp('You have a strong lightning spell. 7 magic per use.')
              elseif deity==3
                  disp('You have a strong earth spell. 7 magic per use.')
              elseif deity==4
                  disp('You have a strong ice spell. 7 magic per use.')
              end
          end
      else
            y_position=y_position;
            x_position=x_position;
      end  %End motion
    
    if current_health<=0
        playing_game=0;
        random_event=24;
    end
    if boss_health<=0;
        playing_game=0;
        random_event=24;
    end
    while x_position<=0  %Boundaries
        x_position=x_position+1;
    end
    while x_position>=20
        x_position=x_position-1;
    end
    while y_position<=0
        y_position=y_position+1;
    end
    while y_position>=20
        y_position=y_position-1;
    end  %finish boundaries
    
    %Now we replot
    clf
    plot(x_home,y_home,'s',x_dungeon,y_dungeon,'p',x_mall,y_mall,'*',x_position,y_position,'d')
    legend('Home','Dungeon','Mall',name)
    axis([0 20 0 20])
    grid on
    
    if x_position==x_home & y_position==y_home  %Begin locations
        current_health=health;
        current_magic=magic;
        disp('Welcome home.  After a good night of rest you feel refreshed.')
    elseif x_position==x_mall & y_position==y_mall
        clc
        shopping=0;
        disp('Welcome to the mall')
        while shopping~=3  %begin shopping
            shopping=input('Where do you want to go?  0- Armory; 1- Weaponry; 2- Item Shop; 3- leave.  ');
            %buy some armor
            if shopping==0  %stores loop
                disp('We have various armors available, which would you like?')
                %warrior and holy knight
                if class==0 | class==27
        armor=input('1- leather armor 10 gold; 2- chain mail 30 gold; 3- plate mail 80 gold   ');
            if armor==1
                if gold >=10
                gold=gold-10
                armorclass=3;
                armor_type=0;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==2
                if gold >=30
                gold=gold-30
                armorclass=5;
                armor_type=1;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==3
                if gold >=80
                gold=gold-80
                armorclass=9;
                armor_type=2;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
                gold=gold
                armmorclass=armorclass;
                armor_type=armor_type;
            end
            %rogue
        elseif class==1
        armor=input('1- cape 5 gold; 2- leather armor 10 gold; 3- light chain mail 40 gold   ');
            if armor==1
                if gold >=5
            gold=gold-5
            armorclass=3;
            armor_type=4;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==2
                if gold >=10
            gold=gold-10
            armorclass=5;
            armor_type=0;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==3
                if gold >=40
            gold=gold-40
            armorclass=7;
            armor_type=5;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
            gold=gold
            armmorclass=armorclass;
                armor_type=armor_type;
            end
            %wizard and priest
        elseif class==2 | class==3
            armor=input('1-robe 5 gold; 2- cloak 15 gold; 3- padded robe 25 gold   ');
             if armor==1
                 if gold >=5
            gold=gold-5
            armorclass=1;
            armor_type=6;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==2
                if gold >=15
            gold=gold-15
            armorclass=2;
            armor_type=7;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif armor==3
                if gold >=25
            gold=gold-25
            armorclass=4;
            armor_type=8;
            else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
            gold=gold
            armmorclass=armorclass;
                armor_type=armor_type;
            end
        end %armors
        %weapons
    elseif shopping==1
        disp('We have weapons for every discipline, take your pick.')
        if class==0 | class==27
            weapon=input('1- long sword 25 gold; 2- mace 30 gold; 3- spear 40 gold   ');
            if weapon==1
                if gold >=25
                gold=gold-25
                weaponclass=6;
                %weapontype is slashing (0), bashing(1) or piercing(2)
                weapontype=0;
                weapon_name=1;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==2
                if gold >=30
                gold=gold-30
                weaponclass=6;
                weapontype=1;
                weapon_name=2;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==3
                if gold >=40
                gold=gold-40
                weaponclass=7;
                weapontype=2;
                weapon_name=3;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
                gold=gold
                weaponclass=weaponclass;
            end
        elseif class==1
            weapon=input('1- knife 10 gold; 2- club 15 gold; 3- rapier 30 gold   ');
            if weapon==1
                if gold >=10
                gold=gold-10
                weaponclass=4;
                weapontype=0;
                weapon_name=4;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==2
                if gold >=15
                gold=gold-15
                weaponclass=4;
                weapontype=1;
                weapon_name=5;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==3   
                if gold >=30
                gold=gold-30
                weaponclass=6;
                weapontype=2;
                weapon_name=6;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
                gold=gold
                weaponclass=weaponclass;
            end
        elseif class==2 | class==3
             weapon=input('1- knife 10 gold; 2- staff 20 gold; 3- dagger 15 gold   ');
            if weapon==1
                if gold >=10
                gold=gold-10
                weaponclass=4;
                weapontype=0;
                weapon_name=4;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==2
                if gold >=20
                    gold=gold-20
                    weaponclass=5;
                    weapontype=1;
                    weapon_name=7;
                    else
                    disp('Sorry, you do not have sufficient funds.')
                end
            elseif weapon==3
                if gold >=15
                    gold=gold-15
                    weaponclass=4;
                    weapontype=2;
                    weapon_name=8;
                else
                    disp('Sorry, you do not have sufficient funds.')
                end
            else
                gold=gold
                weaponclass=weaponclass;
            end%weapons
        end %weapons classes         
            %items
            elseif shopping==2
                disp('Welcome, we have healing herbs and magic tonics, how can we help you?')
                item=input('1- healing herb 10 gold; 2- magic tonic 50 gold  ');
                if item==1
                    if gold >= 10
                        gold=gold-10
                        disp('Thank you for buying a healing herb.')
                        healing_herbs=healing_herbs+1;
                    else
                        disp('Sorry, you do not have sufficient funds.')
                    end
                elseif item==2
                    if gold >=50
                        disp('Thank you for buying a magic tonic.')
                        gold=gold-50
                        magic_tonic=magic_tonic+1;
                    else
                        disp('Sorry, you do not have sufficient funds.')
                    end
                end %item shopping
            end%store choice
    end %store while loop
    elseif x_position==x_dungeon & y_position==y_dungeon
        disp('Welcome to the dungeon')
        location=0;
        random_event=-59;
        clf, clc
        %Here we send it to the dungeon loop
    end %locations all that's open is the location~=0 while
    %random events
    c=ceil(20*rand(1));
    for d=1:c
        r=rand(d);
    end
    orc_health=4;
    %old man
    if random_event==3
        clf
        
        subplot(2,1,1)
        plot(x_home,y_home,'s',x_dungeon,y_dungeon,'p',x_mall,y_mall,'*',x_position,y_position,'d')
        legend('Home','Dungeon','Mall',name)
        axis([0 20 0 20])
        grid on
    
        subplot(2,1,2)
        axis([0 1 0 1])
        axis off
        text(.01, .5, 'You ran into an old man.')
        
        disp('A dirty old man walks up to you starts speaking gibberish.')
        disp('What do you do?')
        old_man_action=input('1- walk away; 2-respond with gibberish; 3- bathe the old man.  ');
        if old_man_action==1
            disp('What a strange experience...')
            experience=experience;
            x_position=x_position+1;
        elseif old_man_action==2
            disp('He looks at you funny and then you share a hearty laugh.')
            experience=experience+ceil(10/level);
        elseif old_man_action==3
            disp('He seems to enjoy it.')
            experience=experience+ceil(20/level);
        else
            disp('What a strange experience...')
            experience=experience;
            x_position=x_position+1;
        end %options for old man
    end %old man event
    %lady
    if random_event==17
        clf
        
        subplot(2,1,1)
        plot(x_home,y_home,'s',x_dungeon,y_dungeon,'p',x_mall,y_mall,'*',x_position,y_position,'d')
        legend('Home','Dungeon','Mall',name)
        axis([0 20 0 20])
        grid on
    
        subplot(2,1,2)
        axis([0 1 0 1])
        axis off
        text(.01, .5, 'You ran into an young lady.')
        
        disp('There is a sexy young lady approaching.  What do you do?')
        lady_action=input('1-walk away; 2- greet her; 3- Try your best pickup line;  ');
        if lady_action==1
            disp('What are you scared of?')
            experience=experience;
            x_position=x_position+1;
        elseif lady_action==2
            disp('She returns your greeting and you have a lovely conversation.')
            experience=experience+ceil(10/level);
        elseif lady_action==3
            pickup_line=input('1- Hey baby, lets integrate! 2- Are you related to my brother?  Do you wanna be?  ');
            acceptance=ceil(2*rand(1));
            if acceptance~=pickup_line
                disp('Slap, she did not care for that one...')
                current_health=current_health-1;
                experience=experience+ceil(25/level);
                if current_health==0
                    current_health=1;
                end
            elseif acceptance==pickup_line
                disp('Nice job it worked.  You the man!')
                current_health=current_health+1;
                experience=experience+ceil(30/level);
                if current_health>health
                    current_health=health;
                end
            else
                disp('Kept your mouth shut, good idea.')
                experience=experience+ceil(5/level);
            end
        else
            disp('What are you scared of?')
            experience=experience;
            x_position=x_position+1;
        end %options for lady
    end %%lady event
    %orc
    if random_event==10  %begin orc event
        clf
        
        subplot(2,1,1)
        plot(x_home,y_home,'s',x_dungeon,y_dungeon,'p',x_mall,y_mall,'*',x_position,y_position,'d')
        legend('Home','Dungeon','Mall',name)
        axis([0 20 0 20])
        grid on
    
        subplot(2,1,2)
        axis([0 1 0 1])
        axis off
        text(.01, .5, 'You ran into an orc.')
        
        disp('You ran into a hostile orc.  What are you gonna do?')
        orc_action=input('1-run away; 2-bribe to go away; 3- FIGHT!!!     ');
        combat_turn=ceil(2*rand(1));  %whose turn
        if orc_action==1  %run away
            disp('You barely made it away.')
            x_position=x_position+1;
        elseif orc_action==2   %Bribe it
            bribe=input('How much are you going to give it?  ');
            if bribe<=gold
                orc_price=round(10*rand(1));
                if bribe>=orc_price
                    disp('The bribe was accepted and it left.')
                    gold=gold-bribe
                else
                    disp('Not enough.  It attacks and you run.')
                    current_health=current_health-3
                end
            else
                disp('You do not have that much gold.  It attacks and you run.')
                current_health=current_health-3
            end            
        elseif orc_action==3  %fight
            orc_health=round(6*rand(1))+round(6*rand(1))+4;
            while orc_health>0 & current_health>0 & combat_turn~=0
                orc_armorclass=4+ceil(4*rand(1));  
                orc_armor_weakness=0;
                orc_element_weakness=1;
                if deity==orc_element_weakness  %bonuses
                    element_bonus=3;
                else
                    element_bonus=0;
                end
                if weapontype==orc_armor_weakness;
                    type_bonus=2;
                else
                    type_bonus=0;
                end
                
                while combat_turn==1
                    disp('Your turn to attack, what are you going to do?')
                    battle_turn=input('1- attack; 2- cure magic; 3- attack magic; 4- use an item; 5- steal; 6- run away.   ');
                        if battle_turn==1  %begin battle turn
                            if orc_armorclass<weaponclass+hitbonus+type_bonus
                                orc_health=orc_health-round(weaponclass*rand(1))-attackbonus-type_bonus;
                                if orc_health<0
                                    orc_health=0;
                                    orcs_action=5;
                                end
                                fprintf('You struck the orc leaving it with %1.0f health.\n',orc_health)
                            else
                                disp('You missed!')
                            end
                        elseif battle_turn==2
                            
                            if cure_option==0
                                disp('You do not know that spell.')
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic>=3
                                disp('You feel much better.')
                                current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                current_magic=current_magic-3;
                                
                                if current_health>health
                                    current_health=health;
                                end
                            
                                fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic<3
                                disp('You do not have enough magic!')
                                
                            elseif cure_option_strong==1
                                cure_strength=input('Which cure will you use? 1- weak cure or 2- strong cure  ');
                                
                                if cure_strength==1 & current_magic>=3
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                    current_magic=current_magic-3;
                                
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==1 & current_magic<3
                                    disp('You do not have enough magic!')
                                    
                                elseif cure_strength==2 & current_magic>=5
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(6*rand(1))+ceil(level*rand(1));
                                    current_magic=current_magic-5;
                                 
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==2 & current_magic<5
                                    disp('You do not have enough magic!')
                                end  %strong cure magic
                            end%cure magic
                            
                        elseif battle_turn==3
                            
                            if attack_option==0
                                disp('You do not know the attack spell!')
                                
                            elseif attack_option==1 & attack_option_strong==0 & current_magic>=4
                                current_magic=current_magic-4;
                                fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %basic cast, finds element
                                if deity==0
                                        disp('You bathe the beast with water!')
                                        orc_health=orc_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        orc_health=orc_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('Cool, you killed it.')
                                         orcs_action=5;
                                     end
                                %short of magic
                          elseif attack_option==1 & attack_option_strong==0 & current_magic<4
                              disp('You are short of magic!')
                          elseif attack_option_strong==1  %strong magic
                                attack_strength=input('How hard will you hit it? 1- weak; 2- strong.  ');
                    %use weak
                                if attack_strength==1 & current_magic>=4
                                    current_magic=current_magic-4;
                                    fprintf('You have %1.0f magic remaining.\n',current_magic)
                             %find element       
                                    if deity==0
                                        disp('You bathe the beast with water!')
                                        orc_health=orc_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        orc_health=orc_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('Cool, its done for.')
                                         orcs_action=5;
                                     end
                                     
                                      %short on magic for weak spell
                                  elseif attack_strength==1 & current_magic<4
                                    disp('You are short of magic!')
                                 %use strong       
                                  elseif attack_strength==2 & current_magic>=7
                                        current_magic=current_magic-7;
                                        fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %finds element
                                        if deity==0
                                            disp('You blast your foe with water!')
                                            orc_health=orc_health-ceil(12*rand(1))-element_bonus-ceil(12*rand(1))-magebonus;
                                        elseif deity==1
                                            disp('You torch your enemy with fire!')
                                            orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                        elseif deity==2
                                              disp('You fry the beast with lightning!')
                                              orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-ceil(12*rand(1))-magebonus;
                                        elseif deity==3
                                              disp('You rock the fiend with a quake!')
                                              orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-ceil(6*rand(1))-ceil(6*rand(1))-magebonus;
                                        elseif deity==4
                                             disp('You freeze the monster with ice!')
                                             orc_health=orc_health-ceil(24*rand(1))-element_bonus-magebonus;
                                       end
                                       fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('It is no more.')
                                         orcs_action=5;
                                     end
                                       %ends elemental attack
                                       %short on magic
                                  elseif attack_strength==2 & current_magic<7
                                        disp('You are short of magic!')
                                  end  %ends strong magic attack
                                  
                              end %ends strong magic option
                              
                          elseif battle_turn==4  %uses an item
                            fprintf('You have %1.0f healing herbs and %1.0f magic tonics.\n',healing_herbs,magic_tonic)
                            item_action=input('What will you use? 1- healing herb; 2- magic tonic?   ');

                            if item_action==1 & healing_herbs>=1                             
                                healing_herbs=healing_herbs-1;
                                disp('You feel better')
                                current_health=current_health+ceil(10*rand(1))+curebonus;
                                
                                if current_health>health
                                  current_health=health;
                                end
                                
                            elseif item_action==1 & healing_herbs<1
                                disp('You do not have any left.')
                                
                            elseif item_action==2 & magic_tonic>=1
                                magic_tonic=magic_tonic-1;
                                current_magic=current_magic+ceil(10*rand(1))+curebonus+magebonus;
                                disp('You feel energized')
                                if current_magic>magic
                                    current_magic=magic;
                                end
                                c=ceil(20*rand(1));
                                  for d=1:c
                                    r=rand(d);
                                    end
                            elseif item_action==2 & magic_tonic<1
                                disp('You do not have any left.')
                            end  %items
                            %steal
                        elseif battle_turn==5
                            steal_chance=ceil(15*rand(1));
                            if steal_option==1
                                if steal_chance<=dexterity/2 & steal_chance>=dexterity/4
                                    disp('You stole a healing herb.')
                                    healing_herbs=healing_herbs+1;
                                elseif steal_chance<=dexterity & steal_chance>dexterity/2
                                    disp('You stole a magic tonic.')
                                    magic_tonic=magic_tonic+1;
                                else
                                    disp('You missed')
                                end
                            else
                                disp('You do not have that skill!')
                            end  %steal
                            
                        elseif battle_turn==6
                           disp('You ran away and lost some gold.')
                           orc_health=0;
                           gold=gold-round(10*rand(1));   
                           if gold<0
                               gold=0;
                           end
                        end %battle_turn
                        if orc_health<=0 | current_health<=0
                            combat_turn=0;
                        else
                            combat_turn=2;
                        end
                    end   
                  while combat_turn==2
                    orcs_action=ceil((10-level)*rand(1));
                    orc_to_hit=ceil(5*rand(1))+3;
                    
                    if orcs_action~=1
                    
                        if orc_to_hit>armorclass
                            orc_attack=round(3*rand(1))+round(3*rand(1));
                            current_health=current_health-orc_attack;
              
                            if current_health>0
                                fprintf('\nThe orc hit you and you now have %1.0f health remaining.\n',current_health)
                            else
                                disp('Sorry, you died')
                            end  %damage done
                            
                        else
                            disp('The orc missed.')
                        end %orc attack
                        
                    elseif orcs_action==1
                        disp('The orc ran away and dropped some gold.')
                        gold=gold+ceil(10*rand(1))
                        orc_health=0;
                    end  %orc actions
                    
                    if orc_health<=0 | current_health<=0  %change turns
                        combat_turn=0;
                    else                        
                        combat_turn=1;
                    end
                    
                end %orc turn
                
                if orc_health<=0 & battle_turn~=6 & orcs_action~=1  %battle's over
                    disp('Congratulations, you won the battle.')
                    experience=experience+ceil((60/level)*rand(1))
                    gold=gold+ceil(20*rand(1))
                end  %you win
                
            end  %ends the while loop
        else
            disp('You barely made it away.')
            x_position=x_position+1;
        end  %ends the combat loop

    end
        %level up
        if experience>=100
            level=level+1;
            experience=experience-100;
            disp('LEVEL UP!')
            if class==0
                health=health+stamina+ceil(6*rand(1))
                magic=magic+ceil(wisdom/3+ceil(2*rand(1)))
                attackbonus=attackbonus+ceil(strength/6);
                hitbonus=hitbonus+ceil(dexterity/6);
                magebonus=magebonus+ceil(intelligence/12);
                curebonus=curebonus+ceil(wisdom/6);
                armorclass=armorclass+2;
                weaponclass=weaponclass+3;
            elseif class==1
                health=health+ceil(stamina/2+ceil(8*rand(1)))
                magic=magic+ceil(wisdom/2+ceil(2*rand(1)))
                attackbonus=attackbonus+ceil(strength/8);
                hitbonus=hitbonus+ceil(dexterity/4);
                magebonus=magebonus+ceil(intelligence/6);
                curebonus=curebonus+ceil(wisdom/12);
                armorclass=armorclass+3;
                weaponclass=weaponclass+2;
           elseif class==2
                health=health+ceil(stamina/4+ceil(4*rand(1)))
                magic=magic+ceil(wisdom+ceil(4*rand(1)))
                attackbonus=attackbonus+ceil(strength/12);
                hitbonus=hitbonus+ceil(dexterity/12);
                magebonus=magebonus+ceil(intelligence/3);
                curebonus=curebonus+ceil(wisdom/4);
                armorclass=armorclass+1;
                weaponclass=weaponclass+1;
            elseif class==3
                health=health+ceil(stamina/3+ceil(5*rand(1)))
                magic=magic+ceil(wisdom+ceil(4*rand(1)))
                attackbonus=attackbonus+ceil(strength/8);
                hitbonus=hitbonus+ceil(dexterity/8);
                magebonus=magebonus+ceil(intelligence/4);
                curebonus=curebonus+ceil(wisdom/3);
                armorclass=armorclass+1;
                weaponclass=weaponclass+2;
            elseif class==27
                health=health+ceil(stamina+ceil(6*rand(1)))
                magic=magic+ceil(wisdom+ceil(4*rand(1)))
                attackbonus=attackbonus+ceil(strength/5);
                hitbonus=hitbonus+ceil(dexterity/5);
                magebonus=magebonus+ceil(intelligence/3);
                curebonus=curebonus+ceil(wisdom/3);
                armorclass=armorclass+3;
                weaponclass=weaponclass+3;
            end  %classes
            %New skills for level up
            
            if level_skill==0
               level_skill=1;
                disp('You have earned a new skill!')
                if cure_option==1 & attack_option==0 & steal_option==0 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong cure; 2- attack magic; 3- stealing ability?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        attack_option=1;
                    elseif new_skill==3
                        steal_option=1;
                    end
                elseif cure_option==1 & attack_option==0 & steal_option==0 & cure_option_strong==1 & attack_option_strong==0
                    new_skill=input('What do you want? 1- attack magic; 2- stealing ability?  ');
                    if new_skill==1
                        attack_option=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                 elseif cure_option==1 & attack_option==1 & steal_option==0 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong cure; 2- stealing ability; 3- strong attack magic?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        steal_option=1;
                    elseif new_skill==3
                        attack_option_strong=1;
                    end  
                  elseif cure_option==1 & attack_option==0 & steal_option==1 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong cure; 2- attack magic?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        attack_option=1;
                    end
                  elseif cure_option==1 & attack_option==1 & steal_option==0 & cure_option_strong==1 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong attack magic; 2- stealing ability?  ');
                    if new_skill==1
                        attack_option_strong=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                  elseif cure_option==1 & attack_option==1 & steal_option==1 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong cure magic; 2- stong attack magic?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        attack_option_strong=1;
                    end
                  elseif cure_option==1 & attack_option==0 & steal_option==0 & cure_option_strong==1 & attack_option_strong==0
                    new_skill=input('What do you want? 1- attack magic; 2- stealing ability?  ');
                    if new_skill==1
                        attack_option=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                  elseif cure_option==1 & attack_option==1 & steal_option==1 & cure_option_strong==1 & attack_option_strong==0
                        disp('You now have strong attack magic.')
                        attack_option_strong=1;
                   elseif cure_option==1 & attack_option==1 & steal_option==0 & cure_option_strong==0 & attack_option_strong==1
                    new_skill=input('What do you want? 1- strong cure; 2- stealing ability?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                 elseif cure_option==1 & attack_option==0 & steal_option==1 & cure_option_strong==1 & attack_option_strong==0
                     disp('You now have a basic attack spell')
                     attack_option=1;
                  elseif cure_option==1 & attack_option==1 & steal_option==0 & cure_option_strong==1 & attack_option_strong==1
                        disp('You now have the ability to steal')
                        steal_option=1;
                  elseif cure_option==0 & attack_option==1 & steal_option==0 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- cure magic; 2- stealing ability; 3- strong attack magic?  ');
                    if new_skill==1
                        cure_option=1;
                    elseif new_skill==2
                        steal_option=1;
                    elseif new_skill==3
                        attack_option_strong=1;
                    end
                  elseif cure_option==0 & attack_option==1 & steal_option==0 & cure_option_strong==0 & attack_option_strong==1
                    new_skill=input('What do you want? 1- cure magic; 2- stealing ability?  ');
                    if new_skill==1
                        cure_option=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                  elseif cure_option==0 & attack_option==1 & steal_option==1 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- cure; 2- strong attack magic?  ');
                    if new_skill==1
                        cure_option=1;
                    elseif new_skill==2
                        attack_option_strong=1;
                    end
                  elseif cure_option==0 & attack_option==1 & steal_option==1 & cure_option_strong==0 & attack_option_strong==1
                    disp('You now have the cure spell.')
                    cure_option=1;
                  elseif cure_option==1 & attack_option==1 & steal_option==1 & cure_option_strong==0 & attack_option_strong==1
                    disp('You now have the strong cure spell')
                    cure_option_strong=1;
                else
                    disp('Sorry, but you already have all available skills.')
                end
            else
                level_skill=0;
            end
           end  %level up
    end  %location town
    c=ceil(20*rand(1));
    for d=1:c
        r=rand(d);
    end
   if location==0;
motion2=input('Use the numberpad to move around the map.   ');
random_event=round(40*rand(1));
    if motion2==1  %begin motion
        x_position2=x_position2-1;
        y_position2=y_position2-1;
    elseif motion2==2
          y_position2=y_position2-1;
          x_position2=x_position2;
    elseif motion2==3
          y_position2=y_position2-1;
          x_position2=x_position2+1;
     elseif motion2==4
          y_position2=y_position2;
          x_position2=x_position2-1;
     elseif motion2==5
          y_position2=y_position2;
          x_position2=x_position2;
     elseif motion2==6
          y_position2=y_position2;
          x_position2=x_position2+1;
     elseif motion2==7
          y_position2=y_position2+1;
          x_position2=x_position2-1;
     elseif motion2==8
          y_position2=y_position2+1;
          x_position2=x_position2;
     elseif motion2==9
          y_position2=y_position2+1;
          x_position2=x_position2+1;
      elseif motion2==0
          y_position2=y_position2;
          x_position2=x_position2;
          random_event=-12;
          clf,clc
          playing_game=0;
          disp('Thanks for playing, good-bye')
      elseif motion2==10
          fprintf('\n%s\nYour current health is %1.0f/%1.0f.\nYour current magic power is %1.0f/%1.0f.\nYour current gold is %1.0f.\nYour current level is %1.0f\nYour current experience is %1.0f.\n',name,current_health,health,current_magic,magic,gold,level,experience)
      elseif motion2==11
          fprintf('You currently have %1.0f gold.\nYou have %1.0f healing herbs.\nYou have %1.0f magic tonics.\n',gold,healing_herbs,magic_tonic)
          if weapon_name==0
              disp('You are fighting bare-handed.')
          elseif weapon_name==1
              disp('You have a long sword.')
          elseif weapon_name==2
              disp('You weapon is a mace.')
          elseif weapon_name==3
              disp('A long spear is your weapon.')
          elseif weapon_name==4
              disp('Your weapon is a knife.')
          elseif weapon_name==5
              disp('You have a club.')
          elseif weapon_name==6
              disp('Your weapon of choice is a rapier.')
          elseif weapon_name==7
              disp('You are carrying a staff.')
          elseif weapon_name==8
              disp('You are holding a dagger.')
          elseif weapon_name==9
              disp('You are weilding the Demonslayer.')
          end
          if armor_type==3
              disp('You have no armor.')
          elseif armor_type==0
              disp('You have leather armor.')
          elseif armor_type==1
              disp('You have chain mail.')
          elseif armor_type==2
              disp('You are donning plate mail.')
          elseif armor_type==4
              disp('You are wearing a cape.')
          elseif armor_type==5
              disp('You have light chain mail.')
          elseif armor_type==6
              disp('You are wearing a robe.')
          elseif armor_type==7
              disp('Your armor is a cloak.')
          elseif armor_type==8
              disp('You have a padded robe.')
          elseif armor_type==9
              disp('You are protected by the Holy Armor.')
          end
           fprintf('\nThe strength of you weapon is %1.0f.\nThe strength of your armor is %1.0f.\n\n',weaponclass,armorclass)
      elseif motion2==12
          fprintf('Your strength is %1.0f.\nYour stamina is %1.0f.\nYour dexterity is %1.0f.\nYour intelligence is %1.0f.\nYour wisdom is %1.0f.\n',strength,stamina,dexterity,intelligence,wisdom)
          if race==0
              disp('You are a human.')
          elseif race==1
              disp('You are an elf.')
          elseif race==2
              disp('You are a dwarf.')
          elseif race==3
              disp('You are a halfling.')
          elseif race==45
              disp('You are a half-elf.')
          end
          if class==0
              disp('You are a warrior.')
          elseif class==1
              disp('You are a rogue.')
          elseif class==2
              disp('You are a wizard.')
          elseif class==3
              disp('You are a priest.')
          elseif class==27
              disp('You are a holy knight.')
          end
          if deity==0
              disp('Your guardian is the serpent Leviathon.')
          elseif deity==1
              disp('Your life is the ash of the Pheonix.')
          elseif deity==2
              disp('The power of Ramuh flows through you.')
          elseif deity==3
              disp('Your strength is that of Titan.')
          elseif deity==4
              disp('Shiva blesses you.')
          end
      elseif motion2==13
          fprintf('Your attack damage bonus is %1.0f.\nYour to hit bonus is %1.0f.\nYour magic attack bonus is %1.0f.\nYour cure magic bonus is %1.0f.\n',attackbonus,hitbonus,magebonus,curebonus)
       elseif motion2==16
          if cure_option==0
             disp('You do not know that spell.')
                               
          elseif cure_option==1 & cure_option_strong==0 & current_magic>=3
                                disp('You feel much better.')
                                current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                current_magic=current_magic-3;
                                
                                if current_health>health
                                    current_health=health;
                                end
                            
                                fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                
          elseif cure_option==1 & cure_option_strong==0 & current_magic<3
                                disp('You do not have enough magic!')
                                
          elseif cure_option_strong==1
                                cure_strength=input('Which cure will you use? 1- weak cure or 2- strong cure  ');
                                
                                if cure_strength==1 & current_magic>=3
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                    current_magic=current_magic-3;
                                
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                 elseif cure_strength==1 & current_magic<3
                                    disp('You do not have enough magic!')
                                    
                                 elseif cure_strength==2 & current_magic>=5
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(6*rand(1))+ceil(level*rand(1));
                                    current_magic=current_magic-5;
                                 
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==2 & current_magic<5
                                   disp('You do not have enough magic!')
                                 end  %strong cure magic
          end%cure magic
      elseif motion2==15
          fprintf('You have %1.0f healing herbs and %1.0f magic tonics.\n',healing_herbs,magic_tonic)
                            item_action=input('What will you use? 1- healing herb; 2- magic tonic?   ');

                            if item_action==1 & healing_herbs>=1                             
                                healing_herbs=healing_herbs-1;
                                disp('You feel better')
                                current_health=current_health+ceil(10*rand(1))+curebonus;
                                
                                if current_health>health
                                  current_health=health;
                                end
                                
                            elseif item_action==1 & healing_herbs<1
                                disp('You do not have any left.')
                                
                            elseif item_action==2 & magic_tonic>=1
                                magic_tonic=magic_tonic-1;
                                current_magic=current_magic+ceil(10*rand(1))+curebonus+magebonus;
                                disp('You feel energized')
                                if current_magic>magic
                                    current_magic=magic;
                                end
                            elseif item_action==2 & magic_tonic<1
                                disp('You do not have any left.')
                            end  %items   
      elseif motion2==14
          if cure_option==1
              disp('You have a basic cure spell. 3 magic per use.')
          end
          if cure_option_strong==1
              disp('You have a strong cure spell. 5 magic per use.')
          end
          if attack_option==1
              if deity==0
                  disp('You have a basic water spell. 4 magic per use.')
              elseif deity==1
                  disp('You have a basic fire spell. 4 magic per use.')
              elseif deity==2
                  disp('You have a basic lightning spell. 4 magic per use.')
              elseif deity==3
                  disp('You have a basic earth spell. 4 magic per use.')
              elseif deity==4
                  disp('You have a basic ice spell. 4 magic per use.')
              end
          end
          if attack_option_strong==1
              if deity==0
                  disp('You have a strong water spell. 7 magic per use.')
              elseif deity==1
                  disp('You have a strong fire spell. 7 magic per use.')
              elseif deity==2
                  disp('You have a strong lightning spell. 7 magic per use.')
              elseif deity==3
                  disp('You have a strong earth spell. 7 magic per use.')
              elseif deity==4
                  disp('You have a strong ice spell. 7 magic per use.')
              end
          end
          c=ceil(20*rand(1));
    for d=1:c
        r=rand(d);
    end
      else
            y_position2=y_position2;
            x_position2=x_position2;
      end  %End motion
       
    if current_health<=0
        playing_game=0;
        random_event=-12;
    end
    if boss_health<=0;
        playing_game=0;
        random_event=-12;
    end
    while x_position2<=0  %Boundaries
        x_position2=x_position2+1;
    end
    while x_position2>=40
        x_position2=x_position2-1;
    end
    while y_position2<=0
        y_position2=y_position2+1;
    end
    while y_position2>=40
        y_position2=y_position2-1;
    end  %finish boundaries
    if random_event==25
        disp('You found a healing fountain.')
        current_health=health;
        current_magic=magic;
    elseif random_event==3 | random_event==17 | random_event==34
        clf
        
        subplot(2,1,1)
        plot(x_position2,y_position2,'d',treasure_location_1x,treasure_location_1y,'p',treasure_location_2x,treasure_location_2y,'p')
        axis([0 40 0 40])
        text(exit_location_x,exit_location_y,'EXIT')
        text(boss_location_x,boss_location_y,'BOSS')
        
        subplot(2,1,2)
        axis([0 1 0 1])
        axis off
        text(.1,.5,'You ran into an item salesman.')
        
        disp('You ran into a traveling item salesman.')
        disp('Welcome, I have healing herbs and magic tonics, how may I help you?')
        item=input('1- healing herb 10 gold; 2- magic tonic 50 gold; 3- nothing;  ');
                if item==1
                    if gold >= 10
                        gold=gold-10
                        disp('Thank you for buying a healing herb.')
                        healing_herbs=healing_herbs+1;
                    else
                        disp('Sorry, you do not have sufficient funds.')
                    end
                elseif item==2
                    if gold >=50
                        disp('Thank you for buying a magic tonic.')
                        gold=gold-50
                        magic_tonic=magic_tonic+1;
                    else
                        disp('Sorry, you do not have sufficient funds.')
                    end
                end %item shopping
    elseif random_event==5 | random_event==39 | random_event==22 | random_event==15 | random_event==27 | random_event==32
        disp('FIGHT')
        clf
        
        subplot(2,1,1)
        plot(x_position2,y_position2,'d',treasure_location_1x,treasure_location_1y,'p',treasure_location_2x,treasure_location_2y,'p')
        axis([0 40 0 40])
        text(exit_location_x,exit_location_y,'EXIT')
        text(boss_location_x,boss_location_y,'BOSS')
        
        subplot(2,1,2)
        axis([0 1 0 1])
        axis off
        text(.1,.5,'You ran into a monster.')
        battle_turn=1;
        monster_generator=ceil(6*rand(1));
        if monster_generator==1 | monster_generator==2 | monster_generator==3
            disp('You ran into an angry marsupial.')
            orc_health=round(8*rand(1))+round(8*rand(1))+8;
            combat_turn=ceil(2*rand(1));
            while orc_health>0 & current_health>0 & combat_turn~=0
                orc_armorclass=6+ceil(6*rand(1));  
                orc_armor_weakness=1;
                orc_element_weakness=4;
                if deity==orc_element_weakness  %bonuses
                    element_bonus=3;
                else
                    element_bonus=0;
                end
                if weapontype==orc_armor_weakness;
                    type_bonus=2;
                else
                    type_bonus=0;
                end
                
                while combat_turn==1
                    disp('Your turn to attack, what are you going to do?')
                    battle_turn=input('1- attack; 2- cure magic; 3- attack magic; 4- use an item; 5- steal; 6- run away.   ');
                        if battle_turn==1  %begin battle turn
                            if orc_armorclass<weaponclass+hitbonus+type_bonus
                                orc_health=orc_health-round(weaponclass*rand(1))-attackbonus-type_bonus;
                                if orc_health<0
                                    orc_health=0;
                                    orcs_action=5;
                                end
                                fprintf('You struck the marsupial leaving it with %1.0f health.\n',orc_health)
                            else
                                disp('You missed!')
                            end
                        elseif battle_turn==2
                            
                            if cure_option==0
                                disp('You do not know that spell.')
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic>=3
                                disp('You feel much better.')
                                current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                current_magic=current_magic-3;
                                
                                if current_health>health
                                    current_health=health;
                                end
                            
                                fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic<3
                                disp('You do not have enough magic!')
                                
                            elseif cure_option_strong==1
                                cure_strength=input('Which cure will you use? 1- weak cure or 2- strong cure  ');
                                
                                if cure_strength==1 & current_magic>=3
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                    current_magic=current_magic-3;
                                
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==1 & current_magic<3
                                    disp('You do not have enough magic!')
                                    
                                elseif cure_strength==2 & current_magic>=5
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(6*rand(1))+ceil(level*rand(1));
                                    current_magic=current_magic-5;
                                 
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==2 & current_magic<5
                                    disp('You do not have enough magic!')
                                end  %strong cure magic
                            end%cure magic
                            
                        elseif battle_turn==3
                            
                            if attack_option==0
                                disp('You do not know the attack spell!')
                                
                            elseif attack_option==1 & attack_option_strong==0 & current_magic>=4
                                current_magic=current_magic-4;
                                fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %basic cast, finds element
                                if deity==0
                                        disp('You bathe the beast with water!')
                                        orc_health=orc_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        orc_health=orc_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('Cool, you killed it.')
                                         orcs_action=5;
                                     end
                                %short of magic
                          elseif attack_option==1 & attack_option_strong==0 & current_magic<4
                              disp('You are short of magic!')
                          elseif attack_option_strong==1  %strong magic
                                attack_strength=input('How hard will you hit it? 1- weak; 2- strong.  ');
                    %use weak
                                if attack_strength==1 & current_magic>=4
                                    current_magic=current_magic-4;
                                    fprintf('You have %1.0f magic remaining.\n',current_magic)
                             %find element       
                                    if deity==0
                                        disp('You bathe the beast with water!')
                                        orc_health=orc_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        orc_health=orc_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('Cool, its done for.')
                                         orcs_action=5;
                                     end
                                      %short on magic for weak spell
                                  elseif attack_strength==1 & current_magic<4
                                    disp('You are short of magic!')
                                 %use strong       
                                  elseif attack_strength==2 & current_magic>=7
                                        current_magic=current_magic-7;
                                        fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %finds element
                                        if deity==0
                                            disp('You blast your foe with water!')
                                            orc_health=orc_health-ceil(12*rand(1))-element_bonus-ceil(12*rand(1))-magebonus;
                                        elseif deity==1
                                            disp('You torch your enemy with fire!')
                                            orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                        elseif deity==2
                                              disp('You fry the beast with lightning!')
                                              orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-ceil(12*rand(1))-magebonus;
                                        elseif deity==3
                                              disp('You rock the fiend with a quake!')
                                              orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-ceil(6*rand(1))-ceil(6*rand(1))-magebonus;
                                        elseif deity==4
                                             disp('You freeze the monster with ice!')
                                             orc_health=orc_health-ceil(24*rand(1))-element_bonus-magebonus;
                                       end
                                       fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('It is no more.')
                                         orcs_action=5;
                                     end
                                       %ends elemental attack
                                       %short on magic
                                  elseif attack_strength==2 & current_magic<7
                                        disp('You are short of magic!')
                                  end  %ends strong magic attack
                                  
                              end %ends strong magic option
                              
                          elseif battle_turn==4  %uses an item
                            fprintf('You have %1.0f healing herbs and %1.0f magic tonics.\n',healing_herbs,magic_tonic)
                            item_action=input('What will you use? 1- healing herb; 2- magic tonic?   ');

                            if item_action==1 & healing_herbs>=1                             
                                healing_herbs=healing_herbs-1;
                                disp('You feel better')
                                current_health=current_health+ceil(10*rand(1))+curebonus;
                                
                                if current_health>health
                                  current_health=health;
                                end
                                
                            elseif item_action==1 & healing_herbs<1
                                disp('You do not have any left.')
                                
                            elseif item_action==2 & magic_tonic>=1
                                magic_tonic=magic_tonic-1;
                                current_magic=current_magic+ceil(10*rand(1))+curebonus+magebonus;
                                disp('You feel energized')
                                if current_magic>magic
                                    current_magic=magic;
                                end
                            elseif item_action==2 & magic_tonic<1
                                disp('You do not have any left.')
                            end  %items
                            %steal
                        elseif battle_turn==5
                            steal_chance=ceil(15*rand(1));
                            if steal_option==1
                                if steal_chance<=dexterity/2 & steal_chance>=dexterity/4
                                    disp('You stole a healing herb.')
                                    healing_herbs=healing_herbs+1;
                                elseif steal_chance<=dexterity & steal_chance>dexterity/2
                                    disp('You stole a magic tonic.')
                                    magic_tonic=magic_tonic+1;
                                else
                                    disp('You missed')
                                end
                            else
                                disp('You do not have that skill!')
                            end  %steal
                            
                        elseif battle_turn==6
                           disp('You ran away and lost some gold.')
                           orc_health=0;
                           gold=gold-round(10*rand(1));   
                           if gold<0
                               gold=0;
                           end
                        end %battle_turn
                        if orc_health<=0 | current_health<=0
                            combat_turn=0;
                        else
                            combat_turn=2;
                        end
                    end   
                  while combat_turn==2
                    orcs_action=ceil((15-level)*rand(1));
                    orc_to_hit=ceil(10*rand(1))+5;
                    
                    if orcs_action~=1
                    
                        if orc_to_hit>armorclass
                            orc_attack=round(5*rand(1))+round(5*rand(1));
                            current_health=current_health-orc_attack;
              
                            if current_health>0
                                fprintf('\nThe marsupial hit you and you now have %1.0f health remaining.\n',current_health)
                            else
                                disp('Sorry, you died')
                            end  %damage done
                            
                        else
                            disp('The marsupial missed.')
                        end %orc attack
                        
                    elseif orcs_action==1
                        disp('The marsupial ran away and dropped some gold.')
                        gold=gold+ceil(25*rand(1))
                        orc_health=0;
                    end  %orc actions
                    
                    if orc_health<=0 | current_health<=0  %change turns
                        combat_turn=0;
                    else                        
                        combat_turn=1;
                    end
                    
                end %orc turn
                
                if orc_health<=0 & battle_turn~=6 & orcs_action~=1  %battle's over
                    disp('Congratulations, you won the battle.')
                    experience=experience+ceil((80/level)*rand(1))
                    gold=gold+ceil(50*rand(1))
                end  %you win
                
            end  %ends the while loop
        
            c=ceil(20*rand(1));
    for d=1:c
        r=rand(d);
    end
        elseif monster_generator==4 | monster_generator==5
            disp('You ran into a hungary hippo.')
            orc_health=round(10*rand(1))+round(10*rand(1))+15;
            combat_turn=ceil(2*rand(1));
            while orc_health>0 & current_health>0 & combat_turn~=0
                orc_armorclass=8+ceil(4*rand(1));  
                orc_armor_weakness=2;
                orc_element_weakness=2;
                if deity==orc_element_weakness  %bonuses
                    element_bonus=3;
                else
                    element_bonus=0;
                end
                if weapontype==orc_armor_weakness;
                    type_bonus=2;
                else
                    type_bonus=0;
                end
                
                while combat_turn==1
                    disp('Your turn to attack, what are you going to do?')
                    battle_turn=input('1- attack; 2- cure magic; 3- attack magic; 4- use an item; 5- steal; 6- run away.   ');
                        if battle_turn==1  %begin battle turn
                            if orc_armorclass<weaponclass+hitbonus+type_bonus
                                orc_health=orc_health-round(weaponclass*rand(1))-attackbonus-type_bonus;
                                if orc_health<0
                                    orc_health=0;
                                    orcs_action=5;
                                end
                                fprintf('You struck the hippo leaving it with %1.0f health.\n',orc_health)
                            else
                                disp('You missed!')
                            end
                        elseif battle_turn==2
                            
                            if cure_option==0
                                disp('You do not know that spell.')
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic>=3
                                disp('You feel much better.')
                                current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                current_magic=current_magic-3;
                                
                                if current_health>health
                                    current_health=health;
                                end
                            
                                fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic<3
                                disp('You do not have enough magic!')
                                
                            elseif cure_option_strong==1
                                cure_strength=input('Which cure will you use? 1- weak cure or 2- strong cure  ');
                                
                                if cure_strength==1 & current_magic>=3
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                    current_magic=current_magic-3;
                                
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==1 & current_magic<3
                                    disp('You do not have enough magic!')
                                    
                                elseif cure_strength==2 & current_magic>=5
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(6*rand(1))+ceil(level*rand(1));
                                    current_magic=current_magic-5;
                                 
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==2 & current_magic<5
                                    disp('You do not have enough magic!')
                                end  %strong cure magic
                            end%cure magic
                            
                        elseif battle_turn==3
                            
                            if attack_option==0
                                disp('You do not know the attack spell!')
                                
                            elseif attack_option==1 & attack_option_strong==0 & current_magic>=4
                                current_magic=current_magic-4;
                                fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %basic cast, finds element
                                if deity==0
                                        disp('You bathe the beast with water!')
                                        orc_health=orc_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        orc_health=orc_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('Cool, you killed it.')
                                         orcs_action=5;
                                     end
                                %short of magic
                          elseif attack_option==1 & attack_option_strong==0 & current_magic<4
                              disp('You are short of magic!')
                          elseif attack_option_strong==1  %strong magic
                                attack_strength=input('How hard will you hit it? 1- weak; 2- strong.  ');
                    %use weak
                                if attack_strength==1 & current_magic>=4
                                    current_magic=current_magic-4;
                                    fprintf('You have %1.0f magic remaining.\n',current_magic)
                             %find element       
                                    if deity==0
                                        disp('You bathe the beast with water!')
                                        orc_health=orc_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        orc_health=orc_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('Cool, its done for.')
                                         orcs_action=5;
                                     end
                                      %short on magic for weak spell
                                  elseif attack_strength==1 & current_magic<4
                                    disp('You are short of magic!')
                                 %use strong       
                                  elseif attack_strength==2 & current_magic>=7
                                        current_magic=current_magic-7;
                                        fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %finds element
                                        if deity==0
                                            disp('You blast your foe with water!')
                                            orc_health=orc_health-ceil(12*rand(1))-element_bonus-ceil(12*rand(1))-magebonus;
                                        elseif deity==1
                                            disp('You torch your enemy with fire!')
                                            orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                        elseif deity==2
                                              disp('You fry the beast with lightning!')
                                              orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-ceil(12*rand(1))-magebonus;
                                        elseif deity==3
                                              disp('You rock the fiend with a quake!')
                                              orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-ceil(6*rand(1))-ceil(6*rand(1))-magebonus;
                                        elseif deity==4
                                             disp('You freeze the monster with ice!')
                                             orc_health=orc_health-ceil(24*rand(1))-element_bonus-magebonus;
                                       end
                                       fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('It is no more.')
                                         orcs_action=5;
                                     end
                                       %ends elemental attack
                                       %short on magic
                                  elseif attack_strength==2 & current_magic<7
                                        disp('You are short of magic!')
                                  end  %ends strong magic attack
                                  
                              end %ends strong magic option
                              
                          elseif battle_turn==4  %uses an item
                            fprintf('You have %1.0f healing herbs and %1.0f magic tonics.\n',healing_herbs,magic_tonic)
                            item_action=input('What will you use? 1- healing herb; 2- magic tonic?   ');

                            if item_action==1 & healing_herbs>=1                             
                                healing_herbs=healing_herbs-1;
                                disp('You feel better')
                                current_health=current_health+ceil(10*rand(1))+curebonus;
                                
                                if current_health>health
                                  current_health=health;
                                end
                                
                            elseif item_action==1 & healing_herbs<1
                                disp('You do not have any left.')
                                
                            elseif item_action==2 & magic_tonic>=1
                                magic_tonic=magic_tonic-1;
                                current_magic=current_magic+ceil(10*rand(1))+curebonus+magebonus;
                                disp('You feel energized')
                                if current_magic>magic
                                    current_magic=magic;
                                end
                            elseif item_action==2 & magic_tonic<1
                                disp('You do not have any left.')
                            end  %items
                            %steal
                        elseif battle_turn==5
                            steal_chance=ceil(15*rand(1));
                            if steal_option==1
                                if steal_chance<=dexterity/2 & steal_chance>=dexterity/4
                                    disp('You stole a healing herb.')
                                    healing_herbs=healing_herbs+1;
                                elseif steal_chance<=dexterity & steal_chance>dexterity/2
                                    disp('You stole a magic tonic.')
                                    magic_tonic=magic_tonic+1;
                                else
                                    disp('You missed')
                                end
                            else
                                disp('You do not have that skill!')
                            end  %steal
                            
                        elseif battle_turn==6
                           disp('You ran away and lost some gold.')
                           orc_health=0;
                           gold=gold-round(10*rand(1));   
                           if gold<0
                               gold=0;
                           end
                        end %battle_turn
                        if orc_health<=0 | current_health<=0
                            combat_turn=0;
                        else
                            combat_turn=2;
                        end
                    end   
                  while combat_turn==2
                    orcs_action=ceil((20-level)*rand(1));
                    orc_to_hit=ceil(15*rand(1))+5;
                    
                    if orcs_action~=1
                    
                        if orc_to_hit>armorclass
                            orc_attack=ceil(7*rand(1))+ceil(7*rand(1))+ceil(7*rand(1));
                            current_health=current_health-orc_attack;
              
                            if current_health>0
                                fprintf('\nThe hippo hit you and you now have %1.0f health remaining.\n',current_health)
                            else
                                disp('Sorry, you died')
                            end  %damage done
                            
                        else
                            disp('The hippo missed.')
                        end %orc attack
                        
                    elseif orcs_action==1
                        disp('The hippo ran away and dropped some gold.')
                        gold=gold+ceil(50*rand(1))
                        orc_health=0;
                    end  %orc actions
                    
                    if orc_health<=0 | current_health<=0  %change turns
                        combat_turn=0;
                    else                        
                        combat_turn=1;
                    end
                    
                end %orc turn
                
                if orc_health<=0 & battle_turn~=6 & orcs_action~=1  %battle's over
                    disp('Congratulations, you won the battle.')
                    experience=experience+ceil((120/level)*rand(1))
                    gold=gold+ceil(150*rand(1))
                end  %you win
                
            end  %ends the while loop
        
            
        elseif monster_generator==6
            disp('You ran into a pissed off panther.')
            orc_health=ceil(15*rand(1))+ceil(15*rand(1))+20;
            combat_turn=ceil(2*rand(1));
            while orc_health>0 & current_health>0 & combat_turn~=0
                orc_armorclass=10+ceil(4*rand(1));  
                orc_armor_weakness=0;
                orc_element_weakness=0;
                if deity==orc_element_weakness  %bonuses
                    element_bonus=3;
                else
                    element_bonus=0;
                end
                if weapontype==orc_armor_weakness;
                    type_bonus=2;
                else
                    type_bonus=0;
                end
                
                while combat_turn==1
                    disp('Your turn to attack, what are you going to do?')
                    battle_turn=input('1- attack; 2- cure magic; 3- attack magic; 4- use an item; 5- steal; 6- run away.   ');
                        if battle_turn==1  %begin battle turn
                            if orc_armorclass<weaponclass+hitbonus+type_bonus
                                orc_health=orc_health-round(weaponclass*rand(1))-attackbonus-type_bonus;
                                if orc_health<0
                                    orc_health=0;
                                    orcs_action=5;
                                end
                                fprintf('You struck the panther leaving it with %1.0f health.\n',orc_health)
                            else
                                disp('You missed!')
                            end
                        elseif battle_turn==2
                            
                            if cure_option==0
                                disp('You do not know that spell.')
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic>=3
                                disp('You feel much better.')
                                current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                current_magic=current_magic-3;
                                
                                if current_health>health
                                    current_health=health;
                                end
                            
                                fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic<3
                                disp('You do not have enough magic!')
                                
                            elseif cure_option_strong==1
                                cure_strength=input('Which cure will you use? 1- weak cure or 2- strong cure  ');
                                
                                if cure_strength==1 & current_magic>=3
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                    current_magic=current_magic-3;
                                
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==1 & current_magic<3
                                    disp('You do not have enough magic!')
                                    
                                elseif cure_strength==2 & current_magic>=5
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(6*rand(1))+ceil(level*rand(1));
                                    current_magic=current_magic-5;
                                 
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==2 & current_magic<5
                                    disp('You do not have enough magic!')
                                end  %strong cure magic
                            end%cure magic
                            
                        elseif battle_turn==3
                            
                            if attack_option==0
                                disp('You do not know the attack spell!')
                                
                            elseif attack_option==1 & attack_option_strong==0 & current_magic>=4
                                current_magic=current_magic-4;
                                fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %basic cast, finds element
                                if deity==0
                                        disp('You bathe the beast with water!')
                                        orc_health=orc_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        orc_health=orc_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('Cool, you killed it.')
                                         orcs_action=5;
                                     end
                                %short of magic
                          elseif attack_option==1 & attack_option_strong==0 & current_magic<4
                              disp('You are short of magic!')
                          elseif attack_option_strong==1  %strong magic
                                attack_strength=input('How hard will you hit it? 1- weak; 2- strong.  ');
                    %use weak
                                if attack_strength==1 & current_magic>=4
                                    current_magic=current_magic-4;
                                    fprintf('You have %1.0f magic remaining.\n',current_magic)
                             %find element       
                                    if deity==0
                                        disp('You bathe the beast with water!')
                                        orc_health=orc_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        orc_health=orc_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('Cool, its done for.')
                                         orcs_action=5;
                                     end
                                      %short on magic for weak spell
                                  elseif attack_strength==1 & current_magic<4
                                    disp('You are short of magic!')
                                 %use strong       
                                  elseif attack_strength==2 & current_magic>=7
                                        current_magic=current_magic-7;
                                        fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %finds element
                                        if deity==0
                                            disp('You blast your foe with water!')
                                            orc_health=orc_health-ceil(12*rand(1))-element_bonus-ceil(12*rand(1))-magebonus;
                                        elseif deity==1
                                            disp('You torch your enemy with fire!')
                                            orc_health=orc_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                        elseif deity==2
                                              disp('You fry the beast with lightning!')
                                              orc_health=orc_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-ceil(12*rand(1))-magebonus;
                                        elseif deity==3
                                              disp('You rock the fiend with a quake!')
                                              orc_health=orc_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-ceil(6*rand(1))-ceil(6*rand(1))-magebonus;
                                        elseif deity==4
                                             disp('You freeze the monster with ice!')
                                             orc_health=orc_health-ceil(24*rand(1))-element_bonus-magebonus;
                                       end
                                       fprintf('\nIt has %1.0f health remaining.\n',orc_health)
                                     if orc_health<=0
                                         disp('It is no more.')
                                         orcs_action=5;
                                     end
                                       %ends elemental attack
                                       %short on magic
                                  elseif attack_strength==2 & current_magic<7
                                        disp('You are short of magic!')
                                  end  %ends strong magic attack
                                  
                              end %ends strong magic option
                              
                          elseif battle_turn==4  %uses an item
                            fprintf('You have %1.0f healing herbs and %1.0f magic tonics.\n',healing_herbs,magic_tonic)
                            item_action=input('What will you use? 1- healing herb; 2- magic tonic?   ');

                            if item_action==1 & healing_herbs>=1                             
                                healing_herbs=healing_herbs-1;
                                disp('You feel better')
                                current_health=current_health+ceil(10*rand(1))+curebonus;
                                
                                if current_health>health
                                  current_health=health;
                                end
                                
                            elseif item_action==1 & healing_herbs<1
                                disp('You do not have any left.')
                                
                            elseif item_action==2 & magic_tonic>=1
                                magic_tonic=magic_tonic-1;
                                current_magic=current_magic+ceil(10*rand(1))+curebonus+magebonus;
                                disp('You feel energized')
                                if current_magic>magic
                                    current_magic=magic;
                                end
                            elseif item_action==2 & magic_tonic<1
                                disp('You do not have any left.')
                            end  %items
                            %steal
                       elseif battle_turn==5
                            steal_chance=ceil(15*rand(1));
                            if steal_option==1
                                if steal_chance<=dexterity/2 & steal_chance>=dexterity/4
                                    disp('You stole a healing herb.')
                                    healing_herbs=healing_herbs+1;
                                elseif steal_chance<=dexterity & steal_chance>dexterity/2
                                    disp('You stole a magic tonic.')
                                    magic_tonic=magic_tonic+1;
                                else
                                    disp('You missed')
                                end
                            else
                                disp('You do not have that skill!')
                            end  %steal
                            
                        elseif battle_turn==6
                           disp('You ran away and lost some gold.')
                           orc_health=0;
                           gold=gold-round(10*rand(1));   
                           if gold<0
                               gold=0;
                           end
                        end %battle_turn
                        if orc_health<=0 | current_health<=0
                            combat_turn=0;
                        else
                            combat_turn=2;
                        end
                    end   
                  while combat_turn==2
                    orcs_action=ceil((50-level)*rand(1));
                    orc_to_hit=ceil(15*rand(1))+10;
                    
                    if orcs_action~=1
                    
                        if orc_to_hit>armorclass
                            orc_attack=ceil(10*rand(1))+ceil(10*rand(1))+ceil(15*rand(1));
                            current_health=current_health-orc_attack;
              
                            if current_health>0
                                fprintf('\nThe panther hit you and you now have %1.0f health remaining.\n',current_health)
                            else
                                disp('Sorry, you died')
                            end  %damage done
                            
                        else
                            disp('The panther missed.')
                        end %orc attack
                        
                    elseif orcs_action==1
                        disp('The panther ran away and dropped some gold.')
                        gold=gold+ceil(25*rand(1))
                        orc_health=0;
                    end  %orc actions
                    
                    if orc_health<=0 | current_health<=0  %change turns
                        combat_turn=0;
                    else                        
                        combat_turn=1;
                    end
                    
                end %orc turn
                
                if orc_health<=0 & battle_turn~=6 & orcs_action~=1  %battle's over
                    disp('Congratulations, you won the battle.')
                    experience=experience+ceil((200/level)*rand(1))
                    gold=gold+ceil(300*rand(1))
                end  %you win
                
            end  %ends the while loop
        end %battles
    elseif random_event==1 | random_event==37
        clf
        
        subplot(2,1,1)
        plot(x_position2,y_position2,'d',treasure_location_1x,treasure_location_1y,'p',treasure_location_2x,treasure_location_2y,'p')
        axis([0 40 0 40])
        text(exit_location_x,exit_location_y,'EXIT')
        text(boss_location_x,boss_location_y,'BOSS')
        
        subplot(2,1,2)
        axis([0 1 0 1])
        axis off
        text(.1,.5,'You found a treasure chest.')
        
        disp('You found a treasure chest.')
        chest_action=input('It may be trapped. 1- try to disarm and open; 2-leave alone;   ');
        if chest_action==1
            trap_chance=round(8*rand(1));
            disarm_chance=dexterity/2;
            if trap_chance<=disarm_chance
                disp('You disarmed the trap!')
                item_chance=ceil(6*rand(1));
                if item_chance==1 | item_chance==2 | item_chance==3
                    disp('You found a healing herb.')
                    healing_herbs=healing_herbs+1;
                elseif item_chance==4 | item_chance==5
                    disp('You found a magic tonic.')
                    magic_tonic=magic_tonic+1;
                elseif item_chance==6
                    disp('You found a healing herb and a magic tonic.')
                    healing_herbs=healing_herbs+1;
                    magic_tonic=magic_tonic+1;
                end
            else
                disp('You failed and the chest exploded!')
                chest_damage=ceil(10*rand(1));
                current_health=current_health-chest_damage;
                if current_health>0
                    fprintf('\nThe chest did %1.0f damage and left you with %1.0f health\n',chest_damage,current_health)
                else
                    fprintf('\n%s died from the explosion.\n',name)
                end
            end
        end
                
    end%random events
    
            
       %level up
        if experience>=100
            level=level+1;
            experience=experience-100;
            disp('LEVEL UP!')
            if class==0
                health=health+stamina+ceil(6*rand(1))
                magic=magic+ceil(wisdom/3+ceil(2*rand(1)))
                attackbonus=attackbonus+ceil(strength/6);
                hitbonus=hitbonus+ceil(dexterity/6);
                magebonus=magebonus+ceil(intelligence/12);
                curebonus=curebonus+ceil(wisdom/6);
                armorclass=armorclass+2;
                weaponclass=weaponclass+3;
            elseif class==1
                health=health+ceil(stamina/2+ceil(8*rand(1)))
                magic=magic+ceil(wisdom/2+ceil(2*rand(1)))
                attackbonus=attackbonus+ceil(strength/8);
                hitbonus=hitbonus+ceil(dexterity/4);
                magebonus=magebonus+ceil(intelligence/6);
                curebonus=curebonus+ceil(wisdom/12);
                armorclass=armorclass+3;
                weaponclass=weaponclass+2;
           elseif class==2
                health=health+ceil(stamina/4+ceil(4*rand(1)))
                magic=magic+ceil(wisdom+ceil(4*rand(1)))
                attackbonus=attackbonus+ceil(strength/12);
                hitbonus=hitbonus+ceil(dexterity/12);
                magebonus=magebonus+ceil(intelligence/3);
                curebonus=curebonus+ceil(wisdom/4);
                armorclass=armorclass+1;
                weaponclass=weaponclass+1;
            elseif class==3
                health=health+ceil(stamina/3+ceil(5*rand(1)))
                magic=magic+ceil(wisdom+ceil(4*rand(1)))
                attackbonus=attackbonus+ceil(strength/8);
                hitbonus=hitbonus+ceil(dexterity/8);
                magebonus=magebonus+ceil(intelligence/4);
                curebonus=curebonus+ceil(wisdom/3);
                armorclass=armorclass+1;
                weaponclass=weaponclass+2;
            elseif class==27
                health=health+ceil(stamina+ceil(6*rand(1)))
                magic=magic+ceil(wisdom+ceil(4*rand(1)))
                attackbonus=attackbonus+ceil(strength/5);
                hitbonus=hitbonus+ceil(dexterity/5);
                magebonus=magebonus+ceil(intelligence/3);
                curebonus=curebonus+ceil(wisdom/3);
                armorclass=armorclass+3;
                weaponclass=weaponclass+3;
            end  %classes
            if level_skill==0
                level_skill=1;
                disp('You have earned a new skill!')
                if cure_option==1 & attack_option==0 & steal_option==0 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong cure; 2- attack magic; 3- stealing ability?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        attack_option=1;
                    elseif new_skill==3
                        steal_option=1;
                    end
                elseif cure_option==1 & attack_option==0 & steal_option==0 & cure_option_strong==1 & attack_option_strong==0
                    new_skill=input('What do you want? 1- attack magic; 2- stealing ability?  ');
                    if new_skill==1
                        attack_option=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                 elseif cure_option==1 & attack_option==1 & steal_option==0 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong cure; 2- stealing ability; 3- strong attack magic?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        steal_option=1;
                    elseif new_skill==3
                        attack_option_strong=1;
                    end  
                  elseif cure_option==1 & attack_option==0 & steal_option==1 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong cure; 2- attack magic?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        attack_option=1;
                    end
                  elseif cure_option==1 & attack_option==1 & steal_option==0 & cure_option_strong==1 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong attack magic; 2- stealing ability?  ');
                    if new_skill==1
                        attack_option_strong=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                  elseif cure_option==1 & attack_option==1 & steal_option==1 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- strong cure magic; 2- stong attack magic?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        attack_option_strong=1;
                    end
                  elseif cure_option==1 & attack_option==0 & steal_option==0 & cure_option_strong==1 & attack_option_strong==0
                    new_skill=input('What do you want? 1- attack magic; 2- stealing ability?  ');
                    if new_skill==1
                        attack_option=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                  elseif cure_option==1 & attack_option==1 & steal_option==1 & cure_option_strong==1 & attack_option_strong==0
                        disp('You now have strong attack magic.')
                        attack_option_strong=1;
                   elseif cure_option==1 & attack_option==1 & steal_option==0 & cure_option_strong==0 & attack_option_strong==1
                    new_skill=input('What do you want? 1- strong cure; 2- stealing ability?  ');
                    if new_skill==1
                        cure_option_strong=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                 elseif cure_option==1 & attack_option==0 & steal_option==1 & cure_option_strong==1 & attack_option_strong==0
                     disp('You now have a basic attack spell')
                     attack_option=1;
                  elseif cure_option==1 & attack_option==1 & steal_option==0 & cure_option_strong==1 & attack_option_strong==1
                        disp('You now have the ability to steal')
                        steal_option=1;
                  elseif cure_option==0 & attack_option==1 & steal_option==0 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- cure magic; 2- stealing ability; 3- strong attack magic?  ');
                    if new_skill==1
                        cure_option=1;
                    elseif new_skill==2
                        steal_option=1;
                    elseif new_skill==3
                        attack_option_strong=1;
                    end
                  elseif cure_option==0 & attack_option==1 & steal_option==0 & cure_option_strong==0 & attack_option_strong==1
                    new_skill=input('What do you want? 1- cure magic; 2- stealing ability?  ');
                    if new_skill==1
                        cure_option=1;
                    elseif new_skill==2
                        steal_option=1;
                    end
                  elseif cure_option==0 & attack_option==1 & steal_option==1 & cure_option_strong==0 & attack_option_strong==0
                    new_skill=input('What do you want? 1- cure; 2- strong attack magic?  ');
                    if new_skill==1
                        cure_option=1;
                    elseif new_skill==2
                        attack_option_strong=1;
                    end
                  elseif cure_option==0 & attack_option==1 & steal_option==1 & cure_option_strong==0 & attack_option_strong==1
                    disp('You now have the cure spell.')
                    cure_option=1;
                  elseif cure_option==1 & attack_option==1 & steal_option==1 & cure_option_strong==0 & attack_option_strong==1
                    disp('You now have the strong cure spell')
                    cure_option_strong=1;
                else
                    disp('Sorry, but you already have all available skills.')
                end
            else
                level_skill=0;
            end
        end  %level up
        clf
      plot(x_position2,y_position2,'d',treasure_location_1x,treasure_location_1y,'p',treasure_location_2x,treasure_location_2y,'p')
      axis([0 40 0 40])
      text(exit_location_x,exit_location_y,'EXIT')
      text(boss_location_x,boss_location_y,'BOSS')
      if x_position2==exit_location_x & y_position2==exit_location_y
          location=1;
      elseif x_position2==boss_location_x & y_position2==boss_location_y
        clf
        
        subplot(2,1,1)
        plot(x_position2,y_position2,'d',treasure_location_1x,treasure_location_1y,'p',treasure_location_2x,treasure_location_2y,'p')
        axis([0 40 0 40])
        text(exit_location_x,exit_location_y,'EXIT')
        text(boss_location_x,boss_location_y,'BOSS')
        
        subplot(2,1,2)
        axis([0 1 0 1])
        axis off
        text(.1,.5,'You ran into Devousius.')
        c=ceil(20*rand(1));
    for d=1:c
        r=rand(d);
    end
          disp('You ran into Devousius.  Prepare to meet your doom.')
            combat_turn=ceil(2*rand(1));
            while boss_health>0 & current_health>0 & combat_turn~=0
                orc_armorclass=25;  
                orc_armor_weakness=0;
                orc_element_weakness=3;
                if deity==orc_element_weakness  %bonuses
                    element_bonus=3;
                else
                    element_bonus=0;
                end
                if weapon_name==9;
                    type_bonus=10;
                else
                    type_bonus=0;
                end
                
                while combat_turn==1
                    disp('Your turn to attack, what are you going to do?')
                    battle_turn=input('1- attack; 2- cure magic; 3- attack magic; 4- use an item; 5- steal; 6- run away.   ');
                        if battle_turn==1  %begin battle turn
                            if orc_armorclass<weaponclass+hitbonus+type_bonus
                                boss_health=boss_health-round(weaponclass*rand(1))-attackbonus-type_bonus;
                                if boss_health<0
                                    boss_health=0;
                                    orcs_action=-5;
                                end
                                fprintf('You struck Devousius leaving him with %1.0f health.\n',boss_health)
                            else
                                disp('You missed!')
                            end
                        elseif battle_turn==2
                            
                            if cure_option==0
                                disp('You do not know that spell.')
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic>=3
                                disp('You feel much better.')
                                current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                current_magic=current_magic-3;
                                
                                if current_health>health
                                    current_health=health;
                                end
                            
                                fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                
                            elseif cure_option==1 & cure_option_strong==0 & current_magic<3
                                disp('You do not have enough magic!')
                                
                            elseif cure_option_strong==1
                                cure_strength=input('Which cure will you use? 1- weak cure or 2- strong cure  ');
                                
                                if cure_strength==1 & current_magic>=3
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(level*rand(1));
                                    current_magic=current_magic-3;
                                
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==1 & current_magic<3
                                    disp('You do not have enough magic!')
                                    
                                elseif cure_strength==2 & current_magic>=5
                                    disp('You feel much better.')
                                    current_health=current_health+round(6*rand(1))+curebonus+ceil(6*rand(1))+ceil(level*rand(1));
                                    current_magic=current_magic-5;
                                 
                                    if current_health>health
                                        current_health=health;
                                    end
                                    
                                    fprintf('\nYour health is now at %1.0f.\nYour magic is now at %1.0f.\n',current_health,current_magic)
                                    
                                elseif cure_strength==2 & current_magic<5
                                    disp('You do not have enough magic!')
                                end  %strong cure magic
                            end%cure magic
                            
                        elseif battle_turn==3
                            
                            if attack_option==0
                                disp('You do not know the attack spell!')
                                
                            elseif attack_option==1 & attack_option_strong==0 & current_magic>=4
                                current_magic=current_magic-4;
                                fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %basic cast, finds element
                                if deity==0
                                        disp('You bathe the beast with water!')
                                        boss_health=boss_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        boss_health=boss_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        boss_health=boss_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        boss_health=boss_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        boss_health=boss_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nDevousius has %1.0f health remaining.\n',boss_health)
                                     if boss_health<=0
                                         disp('You demolished him.')
                                         orcs_action=-5;
                                     end
                                %short of magic
                          elseif attack_option==1 & attack_option_strong==0 & current_magic<4
                              disp('You are short of magic!')
                          elseif attack_option_strong==1  %strong magic
                                attack_strength=input('How hard will you hit it? 1- weak; 2- strong.  ');
                    %use weak
                                if attack_strength==1 & current_magic>=4
                                    current_magic=current_magic-4;
                                    fprintf('You have %1.0f magic remaining.\n',current_magic)
                             %find element       
                                    if deity==0
                                        disp('You bathe the beast with water!')
                                        boss_health=boss_health-ceil(12*rand(1))-element_bonus-magebonus;
                                    elseif deity==1
                                        disp('You burn your enemy with fire!')
                                        boss_health=boss_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                    elseif deity==2
                                        disp('You zap your foe with lightning!')
                                        boss_health=boss_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-magebonus;
                                    elseif deity==3
                                        disp('You smash the monster with rocks!')
                                        boss_health=boss_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-magebonus;
                                    elseif deity==4
                                        disp('You cool the fiend with ice!')
                                        boss_health=boss_health-ceil(10*rand(1))-element_bonus-ceil(2*rand(1))-magebonus;
                                    end
                                    fprintf('\nDevousius has %1.0f health remaining.\n',boss_health)
                                     if boss_health<=0
                                         disp('Cool, he is done for.')
                                         orcs_action=-5;
                                     end
                                      %short on magic for weak spell
                                  elseif attack_strength==1 & current_magic<4
                                    disp('You are short of magic!')
                                 %use strong       
                                  elseif attack_strength==2 & current_magic>=7
                                        current_magic=current_magic-7;
                                        fprintf('You have %1.0f magic remaining.\n',current_magic)
                                %finds element
                                        if deity==0
                                            disp('You blast your foe with water!')
                                            boss_health=boss_health-ceil(12*rand(1))-element_bonus-ceil(12*rand(1))-magebonus;
                                        elseif deity==1
                                            disp('You torch your enemy with fire!')
                                            boss_health=boss_health-ceil(4*rand(1))-element_bonus-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-ceil(4*rand(1))-magebonus;
                                        elseif deity==2
                                              disp('You fry the beast with lightning!')
                                              boss_health=boss_health-ceil(7*rand(1))-element_bonus-ceil(5*rand(1))-ceil(12*rand(1))-magebonus;
                                        elseif deity==3
                                              disp('You rock the fiend with a quake!')
                                              boss_health=boss_health-ceil(6*rand(1))-element_bonus-ceil(6*rand(1))-ceil(6*rand(1))-ceil(6*rand(1))-magebonus;
                                        elseif deity==4
                                             disp('You freeze the monster with ice!')
                                             boss_health=boss_health-ceil(24*rand(1))-element_bonus-magebonus;
                                       end
                                       fprintf('\nDevousius has %1.0f health remaining.\n',boss_health)
                                     if boss_health<=0
                                         disp('He is no more.')
                                         orcs_action=-5;
                                     end
                                       %ends elemental attack
                                       %short on magic
                                  elseif attack_strength==2 & current_magic<7
                                        disp('You are short of magic!')
                                  end  %ends strong magic attack
                                  
                              end %ends strong magic option
                              
                          elseif battle_turn==4  %uses an item
                            fprintf('You have %1.0f healing herbs and %1.0f magic tonics.\n',healing_herbs,magic_tonic)
                            item_action=input('What will you use? 1- healing herb; 2- magic tonic?   ');

                            if item_action==1 & healing_herbs>=1                             
                                healing_herbs=healing_herbs-1;
                                disp('You feel better')
                                current_health=current_health+ceil(10*rand(1))+curebonus;
                                
                                if current_health>health
                                  current_health=health;
                                end
                                
                            elseif item_action==1 & healing_herbs<1
                                disp('You do not have any left.')
                                
                            elseif item_action==2 & magic_tonic>=1
                                magic_tonic=magic_tonic-1;
                                current_magic=current_magic+ceil(10*rand(1))+curebonus+magebonus;
                                disp('You feel energized')
                                if current_magic>magic
                                    current_magic=magic;
                                end
                            elseif item_action==2 & magic_tonic<1
                                disp('You do not have any left.')
                            end  %items
                            %steal
                        elseif battle_turn==5
                            steal_chance=ceil(15*rand(1));
                            if steal_option==1
                                if steal_chance<=dexterity/2 & steal_chance>=dexterity/4
                                    disp('You stole a healing herb.')
                                    healing_herbs=healing_herbs+1;
                                elseif steal_chance<=dexterity & steal_chance>dexterity/2
                                    disp('You stole a magic tonic.')
                                    magic_tonic=magic_tonic+1;
                                else
                                    disp('You missed')
                                end
                            else
                                disp('You do not have that skill!')
                            end  %steal
                            
                        elseif battle_turn==6
                           disp('You ran away and lost some gold.')
                           orc_health=0;
                           gold=gold-round(200*rand(1));   
                           if gold<0
                               gold=0;
                           end
                        end %battle_turn
                        if boss_health<=0 | current_health<=0
                            combat_turn=0;
                        else
                            combat_turn=2;
                        end
                    end   
                  while combat_turn==2
                    orc_to_hit=ceil(25*rand(1))+20;
                    orcs_action=ceil((150-boss_health)*rand(1));
                    if orcs_action>=100
                        disp('Devousius healed himself.')
                        boss_health=ceil(boss_health+40*rand(1)+30);
                        if boss_health>150
                            boss_health=150;
                        end
                    elseif orcs_action>=70 & orcs_action<100
                        current_health=current_health-ceil(10*rand(1))-ceil(10*rand(1))-ceil(10*rand(1))-ceil(10*rand(1))-ceil(10*rand(1))-ceil(10*rand(1))-ceil(10*rand(1))-ceil(10*rand(1))-ceil(10*rand(1))-ceil(10*rand(1));
                        fprintf('\nDevousius hit %s with his breath of death.\n',name)
                         if current_health>0
                                fprintf('\nDevousius hit you and you now have %1.0f health remaining.\n',current_health)
                         else
                                disp('Sorry, you died')
                         end  %damage done
                            
                    else
                        if orc_to_hit>armorclass
                            if holy_armor==1
                                orc_attack=ceil(15*rand(1))+ceil(15*rand(1));
                            else
                                orc_attack=ceil(15*rand(1))+ceil(15*rand(1))+ceil(20*rand(1));
                            end
                            current_health=current_health-orc_attack;
              
                            if current_health>0
                                fprintf('\nDevousius hit you and you now have %1.0f health remaining.\n',current_health)
                            else
                                disp('Sorry, you died')
                            end  %damage done
                            
                        else
                            disp('Devousius missed.')
                        end %orc attack
                    end    
                    if boss_health<=0 | current_health<=0  %change turns
                        combat_turn=0;
                    else                        
                        combat_turn=1;
                    end
                    
                end %orc turn
                
                if boss_health<=0 & battle_turn~=6 %battle's over
                    disp('You destroyed Devousius.')
                    experience=experience+ceil((500/level)*rand(1));
                    gold=gold+ceil(1000*rand(1))
                    if experience>=100
                        level=level+1;
                        experience=experience-100;
                    end
                    experience
                    level
                    clf
                    axis([0 1 0 1])
                    axis off
                    text(.01,.5,'Congratulations, you beat the Dungeon of Devousius!')
                end  %you win
                
            end%end battle
        elseif x_position2==treasure_location_1x & y_position2==treasure_location_1y
            disp('You found a treasure chest.  Inside is the Demonslayer.  Use it to destroy Devousius.')
            weapontype=0;
            weapon_name=9;
            weaponclass=weaponclass+10;
            treasure_location_1x=-10;
            treasure_location_1y=-10;
        elseif x_position2==treasure_location_2x & y_position2==treasure_location_2y
            disp('You found a treasure chest.  Inside is the Holy Armor.  It will protect from Devousius.')
            holy_armor=1;
            armorclass=armorclass+10;
            armor_type=9;
            treasure_location_2x=-10;
            treasure_location_2y=-10;
      end
  end
end%health/playing while
%done